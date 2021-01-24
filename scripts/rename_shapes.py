#!/usr/bin/env python

import argparse
import os
import re
import yaml

parser = argparse.ArgumentParser(description='Rename HH->bbtautau shapes.')
parser.add_argument('--input', required=True, type=str, help="input file")
parser.add_argument('--output', required=True, type=str, help="output file")
parser.add_argument('--config', required=True, type=str, help="configuration file")
args = parser.parse_args()

with open(args.config, 'r') as f:
    cfg = yaml.safe_load(f)

class Filter:
    def __init__(self):
        self.accept = None
        self.reject = None

    def Pass(self, name):
        if self.reject is not None:
            for filter in self.reject:
                if filter.match(name):
                    return False
        if self.accept is not None:
            for filter in self.accept:
                if filter.match(name):
                    return True
            return False
        return True

class RenameRule:
    def __init__(self, pattern, replace):
        self.pattern = re.compile(pattern)
        self.replace = replace

    def Apply(self, name):
        return self.pattern.sub(self.replace, name)

filters = {}
rename_rules = {}

for obj in [ "dir", "hist" ]:
    filter = Filter()
    for sel in [ "accept", "reject" ]:
        f_name = '{}_{}_filters'.format(obj, sel)
        if f_name in cfg:
            f = [ re.compile(s) for s in cfg[f_name] ]
        else:
            f = None
        setattr(filter, sel, f)
    filters[obj] = filter
    if obj + "_rename_rules" in cfg:
        rename_rules[obj] = [ RenameRule(*x) for x in cfg[obj + "_rename_rules"] ]
    else:
        rename_rules[obj] = None


def Pass(name, is_dir):
    obj = "dir" if is_dir else "hist"
    return filters[obj].Pass(name)

def Rename(name, is_dir):
    obj = "dir" if is_dir else "hist"
    if rename_rules[obj] is not None:
        for rule in rename_rules[obj]:
            name = rule.Apply(name)
    return name


import ROOT

f_in = ROOT.TFile(args.input, 'READ')
f_out = ROOT.TFile(args.output, 'RECREATE', '', 209)

for key in sorted(f_in.GetListOfKeys()):
    cl = ROOT.gROOT.GetClass(key.GetClassName())
    if not cl.InheritsFrom("TDirectory"): continue
    if not Pass(key.GetName(), True): continue
    new_name = Rename(key.GetName(), True)
    print('{} -> {}'.format(key.GetName(), new_name))
    dir_in = f_in.Get(key.GetName())
    dir_out = f_out.mkdir(new_name)
    for hist_key in sorted(dir_in.GetListOfKeys()):
        cl = ROOT.gROOT.GetClass(hist_key.GetClassName())
        if not cl.InheritsFrom("TH1"): continue
        if not Pass(hist_key.GetName(), False): continue
        new_hist_name = Rename(hist_key.GetName(), False)
        print('\t{} -> {}'.format(hist_key.GetName(), new_hist_name))
        hist = dir_in.Get(hist_key.GetName()).Clone()
        dir_out.WriteTObject(hist, new_hist_name, "Overwrite")

f_in.Close()
f_out.Close()
