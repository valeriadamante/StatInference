import argparse
import json
import os

parser = argparse.ArgumentParser(description='Rebin shape histograms.')
parser.add_argument('--input', required=True, type=str, help="input directory")
parser.add_argument('--output', required=True, type=str, help="output directory")
parser.add_argument('--bin-edges-json', required=True, type=str, help="json with bin edges")
args = parser.parse_args()

from rebinAndRunLimits import GetLimits

def arrayToStr(a):
    return '[ ' + ', '.join([ str(x) for x in a ]) + ' ]'


with open(args.bin_edges_json, 'r') as f:
    all_bin_edges = json.load(f)

for channel, categories in all_bin_edges.items():
    for category, bin_edges in categories.items():
        file_name = 'hh_{}_{}_13TeV.txt'.format(category, channel)
        input_path = os.path.join(args.input, file_name)
        output_path = os.path.join(args.output, file_name)
        if not os.path.exists(input_path):
            print('{} {}: input not found.'.format(channel, category))
            continue
        if os.path.exists(output_path):
            print('{} {}: output already exists.'.format(channel, category))
            continue
        print('{} {}: rebinning to {} ...'.format(channel, category, arrayToStr(bin_edges)))

        GetLimits(input_path, args.output, bin_edges, 'r', verbose=2, rebin_only=True)
