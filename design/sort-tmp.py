#!/usr/bin/env python3
import os
import glob
import argparse
import numpy as np

parser = argparse.ArgumentParser('Sorting temporary runs.')
parser.add_argument('-d', '--design', type=str, help='Name of the design.')
parser.add_argument('-l', '--model_name', help='Name of the model.')
parser.add_argument('-m', '--mode', choices=['vc', 'cc'],
                    help='Mode of the experiments.')
parser.add_argument('-b', '--best_n', type=int, default=5,
    help='Number of best protocol to be exported.')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

path = f'./out/{args.design}-{args.mode}-{args.model_name}-tmp'

