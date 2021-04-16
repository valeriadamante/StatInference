import bayes_opt
import json
import math
import numpy as np
import os
import queue
from queue import Queue
import random
import re
import ROOT
import shutil
import threading
import time

min_step = 0.001
min_step_digits = -int(math.log10(min_step))
step_int_scale = 10 ** min_step_digits
max_value_int = step_int_scale

def HistToNumpy(hist):
    epsilon = 1e-7
    x = np.zeros(max_value_int)
    if hist.GetNbinsX() != max_value_int:
        raise RuntimeError("Inconsistent number of bins")
    for n in range(max_value_int):
        if abs(hist.GetBinLowEdge(n + 1) - n * min_step) > epsilon \
                or abs(hist.GetBinLowEdge(n + 2) - (n + 1) * min_step) > epsilon:
            raise RuntimeError("Inconsistent binning")
        x[n] = hist.GetBinContent(n + 1)
    return x

def arrayToStr(a):
    return '[ ' + ', '.join([ str(x) for x in a ]) + ' ]'

def ExtractYields(input_shapes, ref_bkgs, nonbkg_regex):
    ref_bkg_yields = {}
    total_bkg_yields = np.zeros(max_value_int)
    for ref_bkg_name in ref_bkgs:
        ref_bkg_yields[ref_bkg_name] = np.zeros(max_value_int)

    input_root = ROOT.TFile.Open(input_shapes)
    hist_names = [ str(key.GetName()) for key in input_root.GetListOfKeys() ]
    nuis_name_regex = re.compile('(.*)_(CMS_.*)(Up|Down)')
    for hist_name in sorted(hist_names):
        if nuis_name_regex.match(hist_name) is not None or nonbkg_regex.match(hist_name) is not None:
            continue
        ref_bkg = None
        for ref_bkg_name, ref_bkg_regex in ref_bkgs.items():
            if ref_bkg_regex.match(hist_name) is not None:
                ref_bkg = ref_bkg_name
        hist = input_root.Get(hist_name)
        bkg_yield = HistToNumpy(hist)
        total_bkg_yields += bkg_yield
        report_str = '{} added to the total bkg contributions'.format(hist_name)
        if ref_bkg is not None:
            ref_bkg_yields[ref_bkg] += bkg_yield
            report_str += ' and {} contributions'.format(ref_bkg)
        print(report_str)
    input_root.Close()
    return ref_bkg_yields, total_bkg_yields


class Binning:
    def __init__(self, edges):
        self.edges = edges
        self.exp_limit = None

    @staticmethod
    def fromRelativeThresholds(rel_thrs, ref_bkg_yields, total_bkg_yields, ref_bkg_thr, total_bkg_thr):
        def ref_bkg_ok(start, stop):
            for bkg_name, bkg_yields in ref_bkg_yields.items():
                if np.sum(bkg_yields[start:stop]) <= ref_bkg_thr:
                    return False
            return True

        def total_bkg_ok(start, stop):
            return np.sum(total_bkg_yields[start:stop]) > total_bkg_thr

        edges = [ max_value_int ]
        for rel_thr in rel_thrs:
            edge_up = edges[-1]
            edge_down = max(edge_up - int(round(rel_thr * step_int_scale)), 0)
            all_ok = False
            while edge_down >= 0:
                all_ok = ref_bkg_ok(edge_down, edge_up) and total_bkg_ok(edge_down, edge_up)
                if all_ok: break
                edge_down -= 1
            if all_ok:
                edges.append(edge_down)
            if edge_down == 0: break
        if len(edges) == 1:
            edges.append(0)
        edges.reverse()
        edges = np.array(edges, dtype=float)
        edges = edges / step_int_scale
        if edges[-1] != 1:
            edges[-1] = 1
        return Binning(edges)

    @staticmethod
    def fromEntry(entry):
        binning = Binning(np.array(entry['bin_edges']))
        binning.exp_limit = float(entry['exp_limit'])
        return binning

    def getRelativeThresholds(self, n_thrs):
        n_inner_edges = len(self.edges) - 2
        if n_thrs < n_inner_edges:
            raise RuntimeError("Invalid number of output thresholds")
        rel_thrs = (np.random.rand(n_thrs) + min_step) / (1 - min_step)
        rel_thrs[0:n_inner_edges] = np.flip(self.edges[2:] - self.edges[1:-1])
        if n_inner_edges < n_thrs:
            rel_thrs[n_inner_edges] = max(self.edges[1], rel_thrs[n_inner_edges])
        return rel_thrs

    def toPoint(self, n_thrs):
        rel_thrs = self.getRelativeThresholds(n_thrs)
        point = {}
        for n in range(n_thrs):
            point['rel_thr_{}'.format(n)] = rel_thrs[n]
        return point

    def isEquivalent(self, other):
        if len(self.edges) != len(other.edges):
            return False
        return np.count_nonzero(self.edges == other.edges) == len(self.edges)
    def isBetter(self, other):
        if other is None or other.exp_limit is None:
            return True
        if self.exp_limit != other.exp_limit:
            return self.exp_limit < other.exp_limit
        return len(self.edges) < len(other.edges)


class BayesianOptimization:
    def __init__(self, max_n_bins, working_area, acq, kappa, xi, input_datacard, poi,
                 ref_bkg_yields, total_bkg_yields, ref_bkg_thr, total_bkg_thr, input_queue_size, random_seed=12345,
                 other_datacards=[]):
        self.max_n_bins = max_n_bins
        self.binnings = []
        self.best_binning = None
        self.working_area = working_area
        self.log_output = os.path.join(working_area, 'results.json')
        self.input_datacard = input_datacard
        self.poi = poi
        self.ref_bkg_yields = ref_bkg_yields
        self.total_bkg_yields = total_bkg_yields
        self.ref_bkg_thr = ref_bkg_thr
        self.total_bkg_thr = total_bkg_thr
        self.other_datacards = other_datacards
        self.best_binning_split = False
        self.input_queue = Queue(input_queue_size)
        self.open_requests = []
        self.binning_lock = threading.Lock()
        self.optimizer_lock = threading.Lock()
        self.print_lock = threading.Lock()

        bounds = {}
        for n in range(max_n_bins):
            bounds['rel_thr_{}'.format(n)] = (min_step, 1.)

        self.bounds_transformer = None
        self.optimizer = bayes_opt.BayesianOptimization(f=None, pbounds=bounds, random_state=random_seed, verbose=1)
        self.utilities = [
            bayes_opt.util.UtilityFunction(kind='ucb', kappa=kappa, xi=xi),# kappa_decay=0.99),
            bayes_opt.util.UtilityFunction(kind='ei', kappa=kappa, xi=xi),
            bayes_opt.util.UtilityFunction(kind='poi', kappa=kappa, xi=xi),
        ]

        if not os.path.isdir(working_area):
            os.mkdir(working_area)
        if os.path.isfile(self.log_output):
            with open(self.log_output, 'r') as f:
                prev_binnings = json.loads('[' + ', '.join(f.readlines()) + ']')
            for entry in prev_binnings:
                binning = Binning.fromEntry(entry)
                if self.findEquivalent(binning) is None:
                    point = binning.toPoint(self.max_n_bins)
                    self.register(point, len(binning.edges), binning.exp_limit)
                    if binning.isBetter(self.best_binning):
                        self.best_binning = binning
                    self.binnings.append(binning)
        if self.best_binning is not None:
            self.print('The best binning from previous iterations: {}, exp_limit = {}' \
                       .format(arrayToStr(self.best_binning.edges), self.best_binning.exp_limit))

        self.bounds_transformer = bayes_opt.SequentialDomainReductionTransformer()
        self.bounds_transformer.initialize(self.optimizer.space)

        suggestions_file = os.path.join(working_area, 'to_try.json')
        self.suggestions = []
        if os.path.isfile(suggestions_file):
            with open(suggestions_file, 'r') as f:
                suggestions = json.load(f)
            for edges in suggestions:
                self.addSuggestion(edges)

    def print(self, msg):
        self.print_lock.acquire()
        print(msg)
        self.print_lock.release()

    def updateBounds(self, n_points):
        if self.bounds_transformer:
            new_bounds = self.bounds_transformer.transform(self.optimizer.space)
            prev_fixed = True
            for n in range(self.max_n_bins):
                key = 'rel_thr_{}'.format(n)
                if prev_fixed:
                    prev_fixed = new_bounds[key][1] - new_bounds[key][0] < 0.1 and n < n_points
                else:
                    new_bounds[key] = np.array([min_step, 1.])
            for key in new_bounds:
                if new_bounds[key][0] >= new_bounds[key][1]:
                    new_bounds[key] = np.array([new_bounds[key][1], new_bounds[key][0]])
                if new_bounds[key][0] < min_step:
                    new_bounds[key][0] = min_step
                if new_bounds[key][1] > 1.:
                    new_bounds[key][1] = 1.
                # if  new_bounds[key][1] - new_bounds[key][0] <= min_step:
                #     new_bounds[key][0] = max(min_step, new_bounds[key][0] - min_step)
                #     new_bounds[key][1] = min(1, new_bounds[key][1] + min_step)
                if np.any(np.isnan(new_bounds[key])):
                    new_bounds[key] = np.array([min_step, 1.])
            self.optimizer.set_bounds(new_bounds)

    def register(self, point, n_edges, exp_limit):
        loss = -(exp_limit*1e6 + n_edges)
        self.optimizer_lock.acquire()
        self.optimizer.register(params=point, target=loss)
        self.updateBounds(n_edges - 2)
        self.optimizer_lock.release()

    def suggest(self, utility_index):
        self.optimizer_lock.acquire()
        if self.suggestions is None or len(self.suggestions) == 0:
            point = self.optimizer.suggest(self.utilities[utility_index])
        else:
            point = self.suggestions[0]
            self.suggestions.remove(point)
        self.optimizer_lock.release()
        return point

    def addSuggestion(self, edges):
        suggested_binning = Binning(np.array(edges))
        rel_thrs = suggested_binning.getRelativeThresholds(self.max_n_bins)
        binning = Binning.fromRelativeThresholds(rel_thrs, self.ref_bkg_yields, self.total_bkg_yields,
                                                 self.ref_bkg_thr, self.total_bkg_thr)
        if self.findEquivalent(binning, False) is None:
            point = binning.toPoint(self.max_n_bins)
            self.suggestions.append(point)

    def findEquivalent(self, binning, lock=True):
        equivalent_binning = None
        if lock:
            self.binning_lock.acquire()
        for prev_binning in self.binnings:
            if prev_binning.isEquivalent(binning):
                equivalent_binning = prev_binning
                break
        if lock:
            self.binning_lock.release()
        return equivalent_binning

    def addNewBinning(self, binning):
        self.binning_lock.acquire()
        if binning.isBetter(self.best_binning):
            self.best_binning = binning
            self.best_binning_split = False
            self.print('New best binning: {}, exp_limit = {}' \
                       .format(arrayToStr(binning.edges), binning.exp_limit))
        self.binnings.append(binning)
        self.binning_lock.release()

    def tryBestBinningSplit(self):
        self.binning_lock.acquire()
        if not self.best_binning_split:
            self.print('Splitting best binning...')
            n_best = len(self.best_binning.edges)
            if n_best - 2 < self.max_n_bins:
                self.optimizer_lock.acquire()
                for k in range(n_best - 1):
                    edges = np.zeros(n_best + 1)
                    edges[0:k+1] = self.best_binning.edges[0:k+1]
                    edges[k+2:] = self.best_binning.edges[k+1:]
                    for delta_edge in [ 0.1, 0.5, 0.9 ]:
                        edges[k+1] = self.best_binning.edges[k] \
                                     + (self.best_binning.edges[k+1] - self.best_binning.edges[k]) * delta_edge
                        self.addSuggestion(edges)
                self.optimizer_lock.release()
            self.best_binning_split = True
        self.binning_lock.release()

    def findOpenRequest(self, binning):
        equivalent_binning = None
        self.binning_lock.acquire()
        for prev_binning in self.open_requests:
            if prev_binning.isEquivalent(binning):
                equivalent_binning = prev_binning
                break
        self.binning_lock.release()
        return equivalent_binning

    def addOpenRequest(self, binning):
        self.binning_lock.acquire()
        self.open_requests.append(binning)
        self.binning_lock.release()

    def clearOpenRequest(self, binning):
        self.binning_lock.acquire()
        to_remove = []
        for prev_binning in self.open_requests:
            if prev_binning.isEquivalent(binning):
                to_remove.append(prev_binning)
        for b in to_remove:
            self.open_requests.remove(b)
        self.binning_lock.release()

    def waitOpenRequestsToFinish(self):
        has_open_requests = True
        while has_open_requests:
            self.binning_lock.acquire()
            has_open_requests = len(self.open_requests) > 0
            self.binning_lock.release()
            time.sleep(1)

    def PointGenerator(self, n_eq_steps):
        n = 0
        utility_index = 0
        while n < n_eq_steps:
            point = self.suggest(utility_index)

            rel_thrs = np.zeros(len(point))
            for k in range(self.max_n_bins):
                rel_thrs[k] = point['rel_thr_{}'.format(k)]
            self.print('rel_thrs: {}'.format(arrayToStr(rel_thrs)))

            binning = Binning.fromRelativeThresholds(rel_thrs, self.ref_bkg_yields, self.total_bkg_yields,
                                                     self.ref_bkg_thr, self.total_bkg_thr)
            equivalent_binning = self.findEquivalent(binning)
            if equivalent_binning is None:
                n = 0
                open_request = self.findOpenRequest(binning)
                if open_request is None:
                    self.print('Next binning to probe: {}'.format(arrayToStr(binning.edges)))
                    self.input_queue.put(binning)
                    self.addOpenRequest(binning)
                    utility_index = 0
                else:
                    self.print('Open request for binning found: {}'.format(arrayToStr(open_request.edges)))
                    time.sleep(1)
                    utility_index = (utility_index + 1) % len(self.utilities)
            else:
                self.print('Equivalent binnig found: {}'.format(arrayToStr(equivalent_binning.edges)))
                self.register(point, len(equivalent_binning.edges), equivalent_binning.exp_limit)
                n += 1
                if n >= 5:
                    self.tryBestBinningSplit()
                if n == n_eq_steps - 1:
                    self.print('Waiting for open requests to finish...')
                    self.waitOpenRequestsToFinish()
        self.input_queue.put(None)

    def JobDispatcher(self):
        while True:
            time.sleep(1)
            for worker_dir in os.listdir(self.working_area):
                worker_dir = os.path.join(self.working_area, worker_dir)
                if not os.path.isdir(worker_dir): continue
                task_file = os.path.join(worker_dir, 'task.txt')
                task_file_tmp = os.path.join(worker_dir, '.task.txt')
                if os.path.isfile(task_file):
                    result_file = os.path.join(worker_dir, 'result.txt')
                    if os.path.isfile(result_file):
                        with open(result_file, 'r') as f:
                            result = json.load(f)
                        binning = Binning.fromEntry(result)
                        self.clearOpenRequest(binning)
                        if self.findEquivalent(binning) is None:
                            point = binning.toPoint(self.max_n_bins)
                            self.addNewBinning(binning)
                            with open(self.log_output, 'a') as f:
                                f.write(json.dumps(result) + '\n')
                            self.register(point, len(binning.edges), binning.exp_limit)

                        os.remove(task_file)
                        os.remove(result_file)
                else:
                    try:
                        binning = self.input_queue.get(True, 1)
                        if binning is None: return
                        task = {
                            'input_datacard': self.input_datacard,
                            'bin_edges': [ x for x in binning.edges ],
                            'poi': self.poi,
                            'other_datacards': self.other_datacards,
                        }
                        with open(task_file_tmp, 'w') as f:
                            json.dump(task, f)
                        shutil.move(task_file_tmp, task_file)
                    except queue.Empty:
                        pass

    def maximize(self, n_eq_steps):

        def generator():
            self.PointGenerator(n_eq_steps)

        def dispatcher():
            self.JobDispatcher()

        generator_thread = threading.Thread(target=generator)
        generator_thread.start()

        dispatcher_thread = threading.Thread(target=dispatcher)
        dispatcher_thread.start()

        generator_thread.join()
        dispatcher_thread.join()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Find optimal binning that minimises the expected limits.')
    parser.add_argument('--input', required=True, type=str, help="input datacard")
    parser.add_argument('--output', required=True, type=str, help="output directory")
    parser.add_argument('--max-n-bins', required=True, type=int, help="maximum number of bins")
    parser.add_argument('--poi', required=False, type=str, default='r', help="parameter of interest")
    parser.add_argument('--verbosity', required=False, type=int, default=1, help="verbosity")
    parser.add_argument('other_datacards', type=str, nargs='*',
                        help="list of other datacards to be combined together with the current target")

    args = parser.parse_args()

    input_datacard = os.path.abspath(args.input)
    input_name = os.path.splitext(os.path.basename(input_datacard))[0]
    input_shapes = os.path.splitext(input_datacard)[0] + '.input.root'

    other_datacards = [ os.path.abspath(p) for p in args.other_datacards ]

    ref_bkgs = {
        'DY': re.compile('^DY.*'),
        'TT': re.compile('^TT$'),
    }

    nonbkg_regex = re.compile('(data_obs|^ggHH.*|^qqHH.*)')

    ref_bkg_yields, total_bkg_yields = ExtractYields(input_shapes, ref_bkgs, nonbkg_regex)
    bo = BayesianOptimization(max_n_bins=args.max_n_bins,
                              working_area=args.output,
                              acq='ucb', kappa=2.57, xi=0.1,
                              input_datacard=input_datacard,
                              poi=args.poi,
                              ref_bkg_yields=ref_bkg_yields,
                              total_bkg_yields=total_bkg_yields,
                              ref_bkg_thr=1e-7, total_bkg_thr=0.18,
                              input_queue_size=2, random_seed=None,
                              other_datacards=other_datacards)
    bo.maximize(20)

    print('Minimization finished.')
    print('Best binning: {}'.format(arrayToStr(bo.best_binning.edges)))
    print('95% CL expected limit: {}'.format(bo.best_binning.exp_limit))
