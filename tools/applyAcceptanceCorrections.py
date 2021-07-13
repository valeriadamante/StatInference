import argparse
import json
import re

parser = argparse.ArgumentParser(description='Apply acceptance corrections.')
parser.add_argument('--input', required=True, type=str, help="input root file")
parser.add_argument('--output', required=True, type=str, help="output root file")
parser.add_argument('--acc-corrs', required=True, type=str, help="acceptance corrections")
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)

processes = set()
unc_sources = set()

f_in = ROOT.TFile(args.input, 'READ')
f_out = ROOT.TFile(args.output, 'RECREATE')

hist_names = sorted([ str(key.GetName()) for key in sorted(f_in.GetListOfKeys()) ])

with open(args.acc_corrs, 'r') as f:
    acc_dict = json.load(f)

corrections = {}

for proc, unc_sources in acc_dict.items():
    for unc_source, unc_scales in unc_sources.items():
        for unc_scale, corr in unc_scales.items():
            hist_name = '{}_{}{}'.format(proc, unc_source, unc_scale)
            corrections[hist_name] = ( proc, corr )

unc_regex = re.compile('.*(Up|Down)$')

central_yeilds = {}

for hist_name in hist_names:
    hist = f_in.Get(hist_name)

    if unc_regex.match(hist_name) is None:
        central_yeilds[hist_name] = hist.Integral()

    new_hist = hist.Clone()
    if hist_name in corrections:
        proc_name, corr = corrections[hist_name]
        old_yield = hist.Integral()
        new_yield = central_yeilds[proc_name] * corr
        if old_yield != 0 and old_yield != new_yield:
            new_hist.Scale(new_yield / old_yield)
            print("{}: {} -> {}".format(hist_name, old_yield, new_yield))

    f_out.WriteTObject(new_hist, hist_name)

f_in.Close()
f_out.Close()
