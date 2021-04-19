import json
import sys

log_file = sys.argv[1]

with open(log_file, 'r') as f:
    binnings = json.loads('[' + ', '.join(f.readlines()) + ']')

def compare_binnings(b1, b2):
    if b2 is None: return True
    if b1['exp_limit'] != b2['exp_limit']: return b1['exp_limit'] < b2['exp_limit']
    if len(b1['bin_edges']) != len(b2['bin_edges']): return len(b1['bin_edges']) < len(b2['bin_edges'])
    for n in reversed(range(len(b1['bin_edges']))):
        if b1['bin_edges'][n] != b2['bin_edges'][n]: return b1['bin_edges'][n] < b2['bin_edges'][n]
    return False

n_bins = []
exp_limits = []

for n in range(21):
    best_binning = None
    for binning in binnings:
        if len(binning['bin_edges']) - 1 == n and compare_binnings(binning, best_binning):
            best_binning = binning
    if best_binning is not None:
        n_bins.append(n)
        exp_limits.append(best_binning['exp_limit'])

def arrayToStr(a):
    return '[ ' + ', '.join([ str(x) for x in a ]) + ' ]'

print('n_bins = {}'.format(arrayToStr(n_bins)))
print('exp_limits = {}'.format(arrayToStr(exp_limits)))

best_binning = None
for binning in binnings:
    if compare_binnings(binning, best_binning):
        best_binning = binning

print("The best binning: {}".format(arrayToStr(best_binning["bin_edges"])))
print("95% CL expected limit: {}".format(best_binning["exp_limit"]))
