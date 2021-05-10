#! /usr/bin/env python3

import argparse
import json
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('path')
parser.add_argument('incoming')
args = parser.parse_args()

incoming_json = json.loads(args.incoming)

try:
  subprocess.check_output(
    [
      'sudo', 'Rscript',
      os.path.join(args.path, 'app', 'main.R'),
      os.path.join(args.path, 'upload'),
      os.path.join(args.path, 'gmt_folder'),
      os.path.join(args.path, 'transform'),
    ]
  )
except subprocess.CalledProcessError as e:
  print(e.output.decode("utf-8"))
  exit(1)

