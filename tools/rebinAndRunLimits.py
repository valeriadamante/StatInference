import datetime
import math
import os
import re
import shutil
import subprocess
import sys
import ROOT


def sh_call(cmd, error_message, verbose=0):
    if verbose > 0:
        print('>> {}'.format(cmd))
    proc = subprocess.Popen(cmd, shell=True, bufsize=1, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = []
    for line in iter(proc.stdout.readline, ""):
        output.append(line)
        if verbose > 1:
            print(line),
    proc.stdout.close()
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError(error_message)
    return output

def ListToVector(x, value_type):
    v = ROOT.std.vector(value_type)(len(x))
    for n in range(len(x)):
        v[n] = x[n]
    return v

def RebinAndFill(new_hist, old_hist):
    epsilon = 1e-7

    def check_range(old_axis, new_axis):
        old_min = old_axis.GetBinLowEdge(1)
        old_max = old_axis.GetBinUpEdge(old_axis.GetNbins())
        new_min = new_axis.GetBinLowEdge(1)
        new_max = new_axis.GetBinUpEdge(new_axis.GetNbins())
        return old_min <= new_min and old_max >= new_max

    def get_new_bin(old_axis, new_axis, bin_id_old):
        old_low_edge = old_axis.GetBinLowEdge(bin_id_old)
        old_up_edge = old_axis.GetBinUpEdge(bin_id_old)
        bin_low_new = new_axis.FindFixBin(old_low_edge)
        bin_up_new = new_axis.FindFixBin(old_up_edge)

        new_up_edge = new_axis.GetBinUpEdge(bin_low_new)
        if not (bin_low_new == bin_up_new or \
                abs(old_up_edge - new_up_edge) <= epsilon * abs(old_up_edge + new_up_edge) * 2):
            raise RuntimeError("Uncompatible bin edges")
        return bin_low_new

    def add_bin_content(bin_old, bin_new):
        cnt_old = old_hist.GetBinContent(bin_old)
        err_old = old_hist.GetBinError(bin_old)
        cnt_new = new_hist.GetBinContent(bin_new)
        err_new = new_hist.GetBinError(bin_new)
        cnt_upd = cnt_new + cnt_old
        err_upd = math.sqrt(err_new ** 2 + err_old ** 2)

        new_hist.SetBinContent(bin_new, cnt_upd)
        new_hist.SetBinError(bin_new, err_upd);

    n_dim = old_hist.GetDimension()
    if new_hist.GetDimension() != n_dim:
        raise RuntimeError("Incompatible number of dimensions")
    if n_dim < 1 or n_dim > 2:
        raise RuntimeError("Unsupported number of dimensions")

    if not check_range(old_hist.GetXaxis(), new_hist.GetXaxis()):
        raise RuntimeError("x ranges are not compatible")

    if n_dim > 1 and not check_range(old_hist.GetYaxis(), new_hist.GetYaxis()):
        raise RuntimeError("y ranges are not compatible")

    for x_bin_old in range(old_hist.GetNbinsX() + 2):
        x_bin_new = get_new_bin(old_hist.GetXaxis(), new_hist.GetXaxis(), x_bin_old)
        if n_dim == 1:
            add_bin_content(x_bin_old, x_bin_new)
        else:
            for y_bin_old in range(old_hist.GetNbinsY() + 2):
                y_bin_new = get_new_bin(old_hist.GetYaxis(), new_hist.GetYaxis(), y_bin_old)
                bin_old = old_hist.GetBin(x_bin_old, y_bin_old)
                bin_new = new_hist.GetBin(x_bin_new, y_bin_new)
                add_bin_content(bin_old, bin_new)

def FixNegativeContributions(histogram):
    correction_factor = 1e-7

    original_integral = histogram.Integral()
    if original_integral <= 0:
        return False

    has_fixed_bins = False
    for n in range(1, histogram.GetNbinsX() + 1):
        if histogram.GetBinContent(n) < 0:
            has_fixed_bins = True
            error = correction_factor - histogram.GetBinContent(n);
            new_error = math.sqrt(error ** 2 + histogram.GetBinError(n) ** 2)
            histogram.SetBinContent(n, correction_factor)
            histogram.SetBinError(n, new_error)

    if has_fixed_bins:
        new_integral = histogram.Integral()
        histogram.Scale(original_integral / new_integral)
    return True

def GetLimits(input_datacard, output_dir, bin_edges, poi, verbose=0):
    input_name = os.path.splitext(os.path.basename(input_datacard))[0]
    input_shapes = os.path.splitext(input_datacard)[0] + '.input.root'
    output_datacard = os.path.join(output_dir, input_name + '.txt')
    output_shapes = os.path.join(output_dir, input_name + '.input.root')

    if verbose > 0:
        print("Preparing datacards and shapes...")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    for out_file in [output_datacard, output_shapes]:
        if os.path.exists(out_file):
            os.remove(out_file)
    shutil.copy(input_datacard, output_datacard)
    input_root = ROOT.TFile.Open(input_shapes)
    output_root = ROOT.TFile(output_shapes, 'RECREATE', '', 209)
    bin_edges_v = ListToVector(bin_edges, 'double')

    processes_to_remove = []
    nuissances_to_remove = []

    hist_names = [ str(key.GetName()) for key in input_root.GetListOfKeys() ]
    name_regex = re.compile('(.*)_(CMS_.*)(Up|Down)')

    for hist_name in sorted(hist_names):
        name_match = name_regex.match(hist_name)
        if name_match is not None:
            process_name = name_match.group(1)
            nuis_name = name_match.group(2)
            is_central = False
        else:
            process_name = hist_name
            is_central = True
        if process_name in processes_to_remove: continue
        hist_orig = input_root.Get(hist_name)
        hist_new = ROOT.TH1F(hist_name, hist_orig.GetTitle(), bin_edges_v.size() - 1, bin_edges_v.data())
        RebinAndFill(hist_new, hist_orig)
        if FixNegativeContributions(hist_new):
            output_root.WriteTObject(hist_new, hist_name)
        else:
            if is_central:
                processes_to_remove.append(process_name)
            else:
                nuissances_to_remove.append('*,{},{}'.format(process_name, nuis_name))

    input_root.Close()
    output_root.Close()

    if len(processes_to_remove):
        proc_str = " ".join(processes_to_remove)
        if verbose > 0:
            print("Removing processes: {}".format(proc_str))
        sh_call('remove_processes.py {} {}'.format(output_datacard, proc_str),
                'Error while running remove_processes.py', verbose)

    if len(nuissances_to_remove):
        nuis_str = " ".join(nuissances_to_remove)
        if verbose > 0:
            print("Removing nuissances: {}".format(nuis_str))
        sh_call('remove_parameters.py {} {}'.format(output_datacard, nuis_str),
                'Error while running remove_parameters.py', verbose)

    if verbose > 0:
        print("Running limits...")
    version = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    law_cmd = 'law run UpperLimits --version {} --datacards {} --pois {} --scan-parameters kl,1,1,1' \
              .format(version, output_datacard, poi)
    output = sh_call(law_cmd, "Error while running UpperLimits", verbose)

    if verbose > 0:
        print("Removing outputs...")
    sh_call(law_cmd + ' --remove-output 2,a ', "Error while removing combine outputs", verbose)
    limit_regex = re.compile('^Expected 50.0%: {} < ([0-9\.]+)'.format(poi))
    for line in reversed(output):
        lim = limit_regex.match(line)
        if lim is not None:
            return float(lim.group(1))

    raise RuntimeError('Limit not found.')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Rebin histogram and run expected limits.')
    parser.add_argument('--input', required=True, type=str, help="input datacard")
    parser.add_argument('--output', required=True, type=str, help="output directory")
    parser.add_argument('--bin-edges', required=True, type=str, help="comma separated bin edges")
    parser.add_argument('--poi', required=False, type=str, default='r', help="parameter of interest")
    parser.add_argument('--verbose', required=False, type=int, default=1, help="verbosity")
    args = parser.parse_args()

    bin_edges = [ float(b) for b in args.bin_edges.split(',') ]
    if args.verbose > 0:
        print('New bin edges: [ {} ]'.format(', '.join([ str(b) for b in bin_edges ])))
    limit = GetLimits(args.input, args.output, bin_edges, args.poi, verbose=args.verbose)

    if args.verbose > 0:
        print('Expected 95% CL limit = {}'.format(limit))
    else:
        print(limit)
