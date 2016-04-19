#!/usr/bin/env python
# encoding: utf-8

import errno
import getopt
import glob
import itertools
import multiprocessing
import os
import sys
import time

import numpy as np

import zscores

probe_set_file_location = 'data/probe-sets'

def make_sure_path_exists(path):
  try:
    os.makedirs(path)
  except OSError as exception:
    if exception.errno != errno.EEXIST:
      raise

def run_zscore(gse, indir, outdir, probe_sets):
  start = time.time()
  path = os.path.join(indir, gse[0], gse[1] + '*.csv')
  probe_set_files = [os.path.join(probe_set_file_location, i, gse[0] + '-' + i.lower().replace('-', '_') + '.csv') for i in probe_sets]
  outpath = os.path.join(outdir, gse[0])
  make_sure_path_exists(outpath)

  gse_file = glob.glob(path)
  if len(gse_file) == 1:
    zscores.script_run(
        gse_file[0],
        gse[2],
        gse[3],
        metadata_feature_rows=gse[4:],
        probe_set_files=probe_set_files,
        outdir=outpath)
  print 'Finished', gse[1], ', took', time.time() - start, 's'

def run_zscore_star(args):
  return run_zscore(args[0], *args[1])

def do_work(indir, outdir, notes_file, probe_sets):
  # Note: columns are
  # 0 GPL
  # 1 GSE
  # 2 Time
  # 3 Censor
  # 4... Optional Metadata Features
  notes = np.genfromtxt(notes_file, dtype=None, delimiter=',')

  pool = multiprocessing.Pool(16)
  remaining_args = [indir, outdir, probe_sets]
  pool.map(run_zscore_star, itertools.izip(notes[1:], itertools.repeat(remaining_args)))

def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'ho:i:n:p:',
    ['help', 'input=', 'output-directory=', 'notes=', 'probe-sets='])
  except getopt.error, msg:
    help_message.usage()

  indir = None
  outdir = '.'
  notes_file = None
  probe_sets = []
  for option, value in opts:
    if option in ('-i', '--input'):
      indir = value
    if option in ('-o', '--output-directory'):
      outdir = value
    if option in ('-n', '--notes'):
      notes_file = value
    if option in ('-p', '--probe-sets'):
      probe_sets = value.split(',')

  if not (indir and notes_file):
    print "Specify input folder and notes file"
    sys.exit(1)
  return indir, outdir, notes_file, probe_sets

def main(argv=None):
  if argv == None:
     argv = sys.argv
  indir, outdir, notes_file, probe_sets = get_options(argv)
  do_work(indir, outdir, notes_file, probe_sets)

if __name__ == '__main__':
  main()


