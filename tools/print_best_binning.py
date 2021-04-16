import json
import sys

log_file = sys.argv[1]

with open(log_file, 'r') as f:
    binnings = json.loads('[' + ', '.join(f.readlines()) + ']')

best_binning = None

def compare_binnings(b1, b2):
    if b2 is None: return True
    if b1['exp_limit'] != b2['exp_limit']: return b1['exp_limit'] < b2['exp_limit']
    if len(b1['bin_edges']) != len(b2['bin_edges']): return len(b1['bin_edges']) < len(b2['bin_edges'])
    for n in reversed(range(len(b1['bin_edges']))):
        if b1['bin_edges'][n] != b2['bin_edges'][n]: return b1['bin_edges'][n] < b2['bin_edges'][n]
    return False

for binning in binnings:
    if compare_binnings(binning, best_binning):
        best_binning = binning

print(best_binning)
