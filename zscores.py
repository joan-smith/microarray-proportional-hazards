#!/usr/bin/env python
# encoding: utf-8
"""
zscores.py

Created by Joan Smith
on 2012-10-13.

Copyright (c) 2015 . All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os
import getopt
import numpy as np
import rpy2.rpy_classic as rpy
from rpy2.robjects.packages import importr
from rpy2.rpy_classic import r as r_old
import rpy2.robjects.numpy2ri
from rpy2.robjects import r
rpy2.robjects.numpy2ri.activate()

help_message = '''This tool allows calculation of Cox Proportional Hazard Z scores and P values for cleaned data files provided
    Usage:
       zscores.py -i <input file> -o <output directory>
    Input file must be provided clean and normalized, though blank data is acceptable.

    Format: a csv file with Unix line endings, each column is a single patient
      rows:
        1. time: first column is an ignored name, all following columns are
             patient survival times
        2. censor: first column is an ignored name, all following columns are
             patient censor values (0 or 1)
        3. patient identifiers: first column is an ignored name, all following
             are patient identifiers
        3...n: first column in a gene/probe name, all following columns are
             patient expression values for that gene/probe
'''

def usage():
  print help_message
  sys.exit(1)

def import_file(name):
  cohort = np.genfromtxt(name, delimiter=',', dtype=None, filling_values='')
  survival_time = [np.float(c) for c in cohort[0][1:] if len(c) > 0]
  survival_censor = [np.int(c) for c in cohort[1][1:] if len(c) > 0]
  gene_names = cohort[:,0][3:]
  patient_value = cohort[3:, 1:]

  for i in patient_value:
    for j in range(len(i)):
      if i[j] == '':
        print 'warning: blank found at row ', i, ' col ', j
        i[j] = np.nan

  patient_value = patient_value.astype(np.float)
  return (survival_time, survival_censor, gene_names, patient_value)

def coxuh(gene_name, expn_value, surv_time, surv_censor):
  rpy.set_default_mode(rpy.NO_CONVERSION)
  r_old.library('survival')

  # remove missing data
  skip_cols = []
  for i in range(len(expn_value)):
    if np.isnan(expn_value[i]):
      skip_cols.append(i)
  expn_value = np.delete(expn_value, skip_cols)
  surv_time = np.delete(surv_time, skip_cols)
  surv_censor = np.delete(surv_censor, skip_cols)

  r.assign('time', surv_time)
  r.assign('censor', surv_censor)
  r.assign('gene', expn_value)
  r('data = data.frame(gene)')
  coxuh_output = r('summary( coxph(formula = Surv(time, censor) ~ gene,' +
    'data = data, model=FALSE, x=FALSE, y=FALSE))')

  coef_ind = list(coxuh_output.names).index('coefficients')
  coeffs = coxuh_output[coef_ind]
  z_ind = list(coeffs.colnames).index('z')
  p_ind = list(coeffs.colnames).index('Pr(>|z|)')
  return {
    'name': gene_name,
    'z': coeffs[z_ind],
    'p': coeffs[p_ind],
    }

def write_file_with_results(input_file_name, results, outfile_location):
  input_file_name_slug = os.path.basename(input_file_name).split('.')[0]
  output_name = os.path.join(outfile_location, input_file_name_slug+'.out.csv')
  with open(output_name, 'w') as outfile:
    outfile.write('Gene/Probe,Z Score, P Value\n')
    for result in results:
      outfile.write(
        result['name'] + ',' +
        '{:g}'.format(result['z']) + ',' +
        '{:g}'.format(result['p']) +
        '\n'
      )

def do_files(files, outdir):
  for f in files:
    survival_time, survival_censor, gene_names, patient_values =  import_file(f)
    results = []
    for i in range(len(patient_values)):
      results.append(coxuh(gene_names[i], patient_values[i] , survival_time , survival_censor))
    write_file_with_results(f, results, outdir)

def main(argv=None):
  if argv is None:
    argv = sys.argv
    try:
      opts, args = getopt.getopt(argv[1:], 'ho:i:v',
      ['help', 'output='])
    except getopt.error, msg:
      usage()

    infile = None
    outdir = '.'
    for option, value in opts:
      if option == '-v':
        verbose = True
      if option in ('-h', '--help'):
        usage()
      if option in ('-i', '--input'):
        infile = value
      if option in ('-o', '--output-directory'):
        outdir = value

    if not infile:
      usage()
    do_files([infile], outdir)

if __name__ == '__main__':
  main()

