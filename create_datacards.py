import os
import sys

if __name__ == "__main__":
  file_dir = os.path.dirname(os.path.abspath(__file__))
  base_dir = os.path.dirname(file_dir)
  file_dir_name = os.path.split(file_dir)[1]
  if base_dir not in sys.path:
    sys.path.append(base_dir)
  __package__ = file_dir_name

from .tools.HarvesterInterface import HarvesterInterface

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='Create datacards.')
  parser.add_argument('--input', required=True, type=str, help="input directory")
  parser.add_argument('--output', required=True, type=str, help="output directory")
  parser.add_argument('--config', required=True, type=str, help="configuration file")
  args = parser.parse_args()

  hi = HarvesterInterface(args.config, args.input)
  hi.createDatacards(args.output)


