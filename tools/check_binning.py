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
from sortedcontainers import SortedSet
import distutils.util

min_step = 0.001
min_step_digits = -int(math.log10(min_step))
step_int_scale = 10 ** min_step_digits
max_value_int = step_int_scale

def HistToNumpy(hist):
    epsilon = 1e-7
    x = np.zeros(max_value_int)
    x_err2 = np.zeros(max_value_int)
    if hist.GetNbinsX() != max_value_int:
        raise RuntimeError("Inconsistent number of bins")
    for n in range(max_value_int):
        if abs(hist.GetBinLowEdge(n + 1) - n * min_step) > epsilon \
                or abs(hist.GetBinLowEdge(n + 2) - (n + 1) * min_step) > epsilon:
            raise RuntimeError("Inconsistent binning")
        x[n] = hist.GetBinContent(n + 1)
        x_err2[n] = hist.GetBinError(n + 1) ** 2
    return x, x_err2

def arrayToStr(a):
    return '[ ' + ', '.join([ str(x) for x in a ]) + ' ]'

class Yields:
    def __init__(self, ref_bkgs, n_bins):
        self.n_bins = n_bins
        self.ref_bkgs = ref_bkgs
        self.yields = {}
        self.processes = SortedSet([ 'total' ] + list(ref_bkgs.keys()))
        self.input_processes = SortedSet()
        self.unc_variations = SortedSet([''])
        self.ref_bkg_thr = 1e-7
        self.total_bkg_thr = 0.03 #0.18
        self.consider_non_central = True
        self.monotonicity_relax_thr = 10.
        self.rel_err_thr = 0.5
        self.rel_err_relax_thr = 1.
        self.rel_unc_variation_thr = 0.2

    def addProcess(self, process, yields, err2, unc_variation):
        names = [ 'total' ]
        for ref_bkg_name, ref_bkg_regex in self.ref_bkgs.items():
            if ref_bkg_regex.match(process) is not None:
                names.append(ref_bkg_name)
        if process not in self.input_processes:
            self.input_processes.add(process)
        for name in names:
            if unc_variation not in self.unc_variations:
                self.unc_variations.add(unc_variation)
            key = (name, unc_variation)
            if key not in self.yields:
                self.yields[key] = [ np.zeros(self.n_bins), np.zeros(self.n_bins) ]
            self.yields[key][0] += yields
            self.yields[key][1] += err2

    def test(self, start, stop, yield_thr):
        total_yield = {}
        print("checking interval [{}, {}] with yield_thr = {}".format(start, stop, yield_thr))
        for process in self.processes:
            is_ref = process != 'total'
            for unc_variation in self.unc_variations:
                is_central = unc_variation == ''
                unc_var_name = 'central' if is_central else unc_variation
                if not is_central and not self.consider_non_central: continue
                key = (process, unc_variation)
                sum = np.sum(self.yields[key][0][start:stop])
                err2 = np.sum(self.yields[key][1][start:stop])
                err = math.sqrt(err2)

                if key not in self.yields:
                    print("Key not found")
                    return False, -1
                if is_ref:
                    thr = self.ref_bkg_thr
                elif is_central:
                    thr = max(self.total_bkg_thr, yield_thr)
                else:
                    thr = self.total_bkg_thr
                if sum < thr:
                    print("{} {}: sum = {} is less than thr = {}. yield_thr = {}".format(process, unc_var_name, sum, thr, yield_thr))
                    return False, -1
                if not is_ref:
                    total_yield[unc_variation] = (sum, math.sqrt(err))
        rel_err = total_yield[''][1] / total_yield[''][0]
        max_delta = 0
        if self.consider_non_central:
            for unc_variation in self.unc_variations:
                max_delta = max(max_delta, abs(total_yield[''][0] - total_yield[unc_variation][0]))
        if rel_err <= self.rel_err_thr or (self.consider_non_central and rel_err <= self.rel_err_relax_thr \
                                           and max_delta / total_yield[''][0] < self.rel_unc_variation_thr):
            print("Interval {} {} accepted".format(start, stop))
            return True, total_yield[''][0]
        print("Interval {} {} failed: rel_err={}".format(start, stop, rel_err))
        return False, -1

    def printSummary(self):
        print('total = {}'.format(' + '.join(self.input_processes)))
        for unc_variation in self.unc_variations:
            unc_name = unc_variation if len(unc_variation) > 0 else 'central'
            print(unc_name)
            for process in self.processes:
                x = 0
                err = 0
                key = (process, unc_variation)
                if key in self.yields:
                    x = np.sum(self.yields[key][0])
                    err = math.sqrt(np.sum(self.yields[key][1]))
                print('\t{}: {:.3f} +/- {:.3f}'.format(process, x, err))

def ExtractYields(input_shapes, ref_bkgs, nonbkg_regex, ignore_variations_regex):
    yields = Yields(ref_bkgs, max_value_int)

    input_root = ROOT.TFile.Open(input_shapes)
    hist_names = [ str(key.GetName()) for key in input_root.GetListOfKeys() ]
    nuis_name_regex = re.compile('(.*)_(CMS_.*(Up|Down))')
    for hist_name in sorted(hist_names):
        if nonbkg_regex.match(hist_name) is not None:
            continue
        nuis_match = nuis_name_regex.match(hist_name)
        if nuis_match is not None:
            process = nuis_match.groups()[0]
            unc_variation = nuis_match.groups()[1]
        else:
            process = hist_name
            unc_variation = ''
        if ignore_variations_regex.match(unc_variation):
            continue

        hist = input_root.Get(hist_name)
        hist_yield, hist_err2 = HistToNumpy(hist)
        yields.addProcess(process, hist_yield, hist_err2, unc_variation)
    input_root.Close()
    return yields

class Binning:
    def __init__(self, edges):
        self.edges = edges
        self.exp_limit = None

    @staticmethod
    def fromRelativeThresholds(rel_thrs, bkg_yields):
        edges = [ max_value_int ]
        prev_yield = -1
        for rel_thr in rel_thrs:
            edge_up = edges[-1]
            edge_down = max(edge_up - int(round(rel_thr * step_int_scale)), 0)
            all_ok = False
            while edge_down >= 0:
                all_ok, new_yield = bkg_yields.test(edge_down, edge_up, prev_yield)
                if all_ok: break
                edge_down -= 1
            if all_ok:
                edges.append(edge_down)
                prev_yield = new_yield
            if edge_down == 0: break
        if len(edges) == 1:
            edges.append(0)
        edges.reverse()
        edges = np.array(edges, dtype=float)
        edges = edges / step_int_scale
        if edges[-1] != 1:
            edges[-1] = 1
        if edges[0] != 0:
            edges[0] = 0
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

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Find optimal binning that minimises the expected limits.')
    parser.add_argument('--input', required=True, type=str, help="input datacard")
    parser.add_argument('--edges', required=True, type=str, help="bin edges")
    parser.add_argument('--params', required=False, type=str, default=None,
                        help="algorithm parameters in format param1=value1,param2=value2 ...")

    args = parser.parse_args()

    input_datacard = os.path.abspath(args.input)
    input_name = os.path.splitext(os.path.basename(input_datacard))[0]
    input_shapes = os.path.splitext(input_datacard)[0] + '.input.root'

    param_dict = {}
    if args.params is not None:
        params = args.params.split(',')
        for param in params:
            p = param.split('=')
            if len(p) != 2:
                raise RuntimeError('invalid parameter definition "{}"'.format(param))
            param_dict[p[0]] = p[1]

    ref_bkgs = {
        'DY': re.compile('^DY$'),
        'TT': re.compile('^TT$'),
    }

    nonbkg_regex = re.compile('(data_obs|^ggHH.*|^qqHH.*|^DY_[0-2]b.*)')
    ignore_unc_variations = re.compile('(CMS_bbtt_201[6-8]_DYSFunc[0-9]+|CMS_bbtt_.*_QCDshape)(Up|Down)')
    print("Extracting yields for background processes...")
    bkg_yields = ExtractYields(input_shapes, ref_bkgs, nonbkg_regex, ignore_unc_variations)
    #bkg_yields.printSummary()

    for name,value in param_dict.items():
        if not hasattr(bkg_yields, name):
            raise RuntimeError('Unknown parameter "{}"'.format(name))
        def_value = getattr(bkg_yields, name)
        def_value_type = type(def_value)
        if def_value_type == bool:
            new_value = bool(distutils.util.strtobool(value))
        else:
            new_value = def_value_type(value)
        print("Setting {} = {}".format(name, new_value))
        setattr(bkg_yields, name, new_value)

    max_n_bins = 20
    edges = json.loads('[{}]'.format(args.edges))
    suggested_binning = Binning(np.array(edges))
    rel_thrs = suggested_binning.getRelativeThresholds(max_n_bins)
    binning = Binning.fromRelativeThresholds(rel_thrs, bkg_yields)

    print("Resulting binning: {}".format(arrayToStr(binning.edges)))
