import argparse
import json
import re

parser = argparse.ArgumentParser(description='Extract integrals into json shape variations.')
parser.add_argument('--input', required=True, type=str, help="input root file")
parser.add_argument('--output', required=True, type=str, help="output json file")
parser.add_argument('--bin', required=True, type=str, help="bin")
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)

yields = {}

f_in = ROOT.TFile(args.input, 'READ')
bin_dir = f_in.Get(args.bin)

unc_regex = re.compile('.*(Up|Down)$')

hist_names = [ str(key.GetName()) for key in sorted(bin_dir.GetListOfKeys()) ]
for hist_name in hist_names:
    unc_match = unc_regex.match(hist_name)
    if unc_match is None:
        if hist_name in yields:
            raise RuntimeError('Duplicated process {}'.format(hist_name))
        yields[hist_name] = {}
        hist = bin_dir.Get(hist_name)
        yields[hist_name]['Central'] = hist.Integral()

for hist_name in hist_names:
    unc_match = unc_regex.match(hist_name)
    if unc_match is not None:
        best_proc = None
        best_match = None
        for proc in yields.keys():
            unc_src_match = re.match('^{}_(.*)(Up|Down)$'.format(proc), hist_name)
            if unc_src_match is not None:
                if best_proc is not None and len(best_proc) == len(proc):
                    raise RuntimeError('Two processes match {}: {} and {}'.format(hist_name, best_proc, proc))
                if best_proc is None or len(best_proc) < len(proc):
                    best_proc = proc
                    best_match = unc_src_match
        if best_proc is None:
            raise RuntimeError('Process not found for {}'.format(hist_name))
        unc_source = best_match.groups()[0]
        unc_scale = best_match.groups()[1]
        if unc_source not in yields[best_proc]:
            yields[best_proc][unc_source] = {}
        if unc_scale in yields[best_proc][unc_source]:
            raise RuntimeError('Duplicated unc_scale={} for {} {}'.format(unc_scale, unc_source, best_proc))
        hist = bin_dir.Get(hist_name)
        yields[best_proc][unc_source][unc_scale] = hist.Integral()

f_in.Close()


with open(args.output, 'w') as f:
    json.dump(yields, f, sort_keys=True, indent=4)
