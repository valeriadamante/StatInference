import argparse
import json
import re

parser = argparse.ArgumentParser(description='Extract integrals into json shape variations.')
parser.add_argument('input', nargs='+', type=str, help="input root files")
#parser.add_argument('--bin', required=False, default=None, type=str, help="bin")
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)

processes = set()
unc_sources = set()

for file_name in args.input:
    f_in = ROOT.TFile(file_name, 'READ')
    #if args.bin is None:
    bin_dir = f_in
    # else:
    #     bin_dir = f_in.Get(args.bin)

    unc_regex = re.compile('.*(Up|Down)$')

    hist_names = [ str(key.GetName()) for key in sorted(bin_dir.GetListOfKeys()) ]

    f_in.Close()

    for hist_name in hist_names:
        unc_match = unc_regex.match(hist_name)
        if unc_match is None:
            processes.add(hist_name)


    for hist_name in hist_names:
        unc_match = unc_regex.match(hist_name)
        if unc_match is not None:
            best_proc = None
            best_match = None
            for proc in processes:
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
            unc_sources.add(unc_source)

print('processes = {')
for process in sorted(processes):
    print('\t "{0}" : "{0}",'.format(process))
print('}')

print('unc_sources = {')
for unc_source in sorted(unc_sources):
    print('\t "{0}" : "{0}",'.format(unc_source))
print('}')
