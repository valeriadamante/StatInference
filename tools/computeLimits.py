import argparse
import datetime
import json
import os
import re
import subprocess

parser = argparse.ArgumentParser(description='Rebin shape histograms.')
parser.add_argument('--input', required=True, type=str, help="input directory")
parser.add_argument('--channel', required=True, type=str, help="channel_year")
args = parser.parse_args()


def arrayToStr(a):
    return '[ ' + ', '.join([ str(x) for x in a ]) + ' ]'

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

def extractLimits(output, poi):
    limit_regex = re.compile('^Expected 50.0%: {} < ([0-9\.]+)'.format(poi))
    for line in reversed(output):
        lim = limit_regex.match(line)
        if lim is not None:
            return float(lim.group(1))
    raise RuntimeError('Limit not found.')



output_json = os.path.join(args.input, '{}_limits.json'.format(args.channel))
if os.path.exists(output_json):
    with open(output_json, 'r') as f:
        limits = json.load(f)
else:
    limits = {}

#categories = [ 'res2b', 'res1b', 'classVBF', 'classGGF', 'classttH', 'classTT', 'classDY', 'boosted' ]
categories = [ 'res2b', 'res1b', 'boosted', 'classVBF', 'classGGF', 'classttH', 'classTT', 'classDY' ]
pois = [ 'r', 'r_qqhh' ]

for n in range(len(categories)):
    cat_cmb = '_'.join(categories[:n+1])

    if cat_cmb in limits:
        missing_pois = []
        for poi in pois:
            if poi not in limits[cat_cmb]:
                missing_pois.append(poi)
    else:
        limits[cat_cmb] = {}
        missing_pois = pois
    if len(missing_pois) > 0:
        version = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        datacards = []
        has_all_datacards = True
        for cat in categories[:n+1]:
            datacard = os.path.join(args.input, 'hh_{}_{}_13TeV.txt'.format(cat, args.channel))
            if not os.path.exists(datacard):
                has_all_datacards = False
                break
            datacards.append(datacard)
        if not has_all_datacards:
            break
        datacards_str = ','.join(datacards)
        for poi in missing_pois:
            law_cmd = 'law run UpperLimits --version {} --hh-model {} --datacards {} --pois {} --scan-parameters {}' \
                      .format(version, 'hh_model.model_default', datacards_str, poi, 'kl,1,1,1')
            output = sh_call(law_cmd, "Error while running UpperLimits", 2)
            limits[cat_cmb][poi] =  extractLimits(output, poi)

        sh_call(law_cmd + ' --remove-output 2,a ', "Error while removing combine outputs", 2)

    print("{} limits: r = {}, r_qqhh = {}".format(cat_cmb, limits[cat_cmb]['r'], limits[cat_cmb]['r_qqhh']))
    with open(output_json, 'w') as f:
        json.dump(limits, f, indent=4)
