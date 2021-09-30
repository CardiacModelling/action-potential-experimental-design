#!/usr/bin/env python3
import os
import glob
import re
import argparse
import numpy as np
from subprocess import call

parser = argparse.ArgumentParser('Sorting temporary runs.')
parser.add_argument('-d', '--design', type=str, help='Name of the design.')
parser.add_argument('-l', '--model_name', help='Name of the model.')
parser.add_argument('-m', '--mode', choices=['vc', 'cc'],
                    help='Mode of the experiments.')
parser.add_argument('-b', '--best_n', type=int, default=5,
    help='Number of best protocol to be exported.')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

run = 0

path = f'./out/{args.design}-{args.mode}-{args.model_name}-tmp'
new_path = f'./out/{args.design}-{args.mode}-{args.model_name}'
if not os.path.isdir(new_path):
    os.makedirs(new_path)

protocols = glob.glob(f'{path}/opt-prt-run{run}-re[0-9]-[0-9].txt')
scores = glob.glob(f'{path}/opt-prt-score-run{run}-re[0-9]-[0-9].txt')
assert(type(protocols) is list)
assert(len(protocols) == len(scores))

score_values = []
for protocol in protocols:
    i, j = re.findall(f'{path}/opt-prt-run{run}-re(\d)-(\d).txt', protocol)[0]
    s = f'{path}/opt-prt-score-run{run}-re{i}-{j}.txt'
    assert(s in scores)
    score = np.loadtxt(s)
    score_values.append(score)

call(['cp',
      f'{path}/opt-prt-run{run}-re{i}.out',
      f'{new_path}/opt-prt-run{run}.out'])

call(['cp',
      f'{path}/opt-prt-run{run}-re{i}-parameter_evaluate.txt',
      f'{new_path}/opt-prt-run{run}-parameter_evaluate.txt'])

order = np.argsort(score_values)  # (use [::-1] for LL)

for i, o in enumerate(order[:args.best_n]):
    protocol = protocols[o]
    score = score_values[o]
    call(['cp',
          protocol,
          f'{new_path}/opt-prt-run{run}-rank{o}.txt'])
    np.savetxt(f'{new_path}/opt-prt-score-run{run}-rank{o}.txt', [score])
