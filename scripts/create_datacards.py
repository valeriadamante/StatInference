#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description='Create HH->bbtautau datacards.')
parser.add_argument('--input', required=True, type=str, help="input directory")
parser.add_argument('--output', required=True, type=str, help="output directory")
parser.add_argument('--config', required=True, type=str, help="configuration file")
args = parser.parse_args()

from tools.HarvesterInterface import HarvesterInterface

hi = HarvesterInterface(args.config, args.input)

for year, channel, category in hi.YCC():
    hi.AddObservations(year, channel, category)
hi.cb.PrintAll()
