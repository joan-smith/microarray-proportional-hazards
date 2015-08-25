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

import pdb
import sys
import re
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
       ./zscores.py -i <input file> -o <output directory>
          -m <feature_name>,<feature_name>|all
          --interactive

    Input file must be provided clean and normalized, though blank data is acceptable.

    Format: a csv file with Unix line endings, each column is a single patient
    In interactive mode, all fields except ID_REF are optional
      rows:
        1. time: first column is an ignored name, all following columns are
             patient survival times
        2. censor: first column is an ignored name, all following columns are
             patient censor values (0 or 1)
        3..n features: (optional) features for multivariate analysis. Features are selected by their name,
             in the first column.
        n + 1. patient identifiers: first column is "patient", all following
             are patient identifiers
        n + 2 ...m : first column in a gene/probe name, all following columns are
             patient expression values for that gene/probe
'''

def usage():
  print help_message
  sys.exit(1)

def safe_string(string):
  safe_string = "".join([c for c in string if re.match(r'\w', c)])
  return safe_string

def print_row_selection(row_titles):
  for i, title in enumerate(row_titles):
    print str(i+1) + '. ' + title

def row_selection(row_selection, row_titles, exit_on_invalid = True):
  try:
    row_selection_int = int(row_selection) - 1
    if row_selection_int >= len(row_titles):
      print 'Error: Invalid Selection'
      if exit_on_invalid:
        exit(1)
      return None
    return row_selection_int
  except(ValueError):
    print 'Error: Invalid Selection'
    if exit_on_invalid:
      exit(1)
    return None

def import_file_interactive(name):
  with open(name) as f:
    row_titles = []
    for line in f.readlines():
      row_title = line.split(',')[0]
      if row_title == 'patient':
        break
      else:
        row_titles.append(row_title)

    print_row_selection(row_titles)
    time_row_number_selection = raw_input('Enter number for row that has survival time: ')
    time_row_number = row_selection(time_row_number_selection, row_titles)

    print_row_selection(row_titles)
    censor_row_number_selection = raw_input('Enter number for row that has censor: ')
    censor_row_number = row_selection(censor_row_number_selection, row_titles)

    additional_variables_rows = []
    while True:
      print_row_selection(row_titles)
      selection = raw_input('Enter number for row with additional value, or "END" to finish selection: ')
      if selection == 'END':
        break
      else:
        additional_variables_rows.append(row_selection(selection, row_titles, exit_on_invalid=False))
    additional_variables_rows = [x for x in additional_variables_rows if x != None]
    return import_file(name, time_row_number, censor_row_number, features_rows=additional_variables_rows)

def import_file(name, time_row=0, censor_row=1, features_rows=None):
  cohort = np.genfromtxt(name, delimiter=',', dtype=None, filling_values='')
  survival_time = [np.float(c) for c in cohort[time_row][1:] if len(c) > 0]
  survival_censor = [np.int(c) for c in cohort[censor_row][1:] if len(c) > 0]
  row_headers = list(cohort[:,0])

  if 'patient' in row_headers:
    patient_row_idx = row_headers.index('patient')
  else:
    print 'Error: \'patient\' header not found'
    usage()

  if features_rows == None:
    feature_names = [safe_string(name) for name in cohort[2:patient_row_idx,0]]
    features = cohort[2:patient_row_idx,1:]
  else:
    feature_names = [safe_string(name) for name in cohort[features_rows, 0]]
    features = cohort[features_rows,1:]

  for i, feature in enumerate(features):
    for j in range(len(feature)):
      if feature[j] == '':
        feature[j] = np.nan

  gene_names = row_headers[patient_row_idx+1:]
  patient_value = cohort[patient_row_idx+1:, 1:]
  for i, row in enumerate(patient_value):
    for j in range(len(row)):
      if row[j] == '' or row[j] == np.nan:
        row[j] = np.nan

  patient_value = patient_value.astype(np.float)
  return (survival_time, survival_censor, gene_names, patient_value, feature_names, features)

def coxuh(gene_name, expn_value, surv_time, surv_censor, feature_names, features):
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
  if len(feature_names) >= 1:
    features = np.delete(features, skip_cols, 1)

  if len(expn_value) < 10:
    print 'warning: not enough samples for row ' + gene_name + ', skipping.'
    return {}

  r.assign('time', surv_time)
  r.assign('censor', surv_censor)
  for idx, feature_name in enumerate(feature_names):
    float_features = features[idx].astype(np.float)
    r.assign(feature_name, float_features)
  formula_string = ''
  if len(feature_names) >= 1:
    formula_string = 'gene + ' + ' + '.join(feature_names)
    data_frame_string = 'gene, '+ ', '.join(feature_names)
  else:
    formula_string = 'gene'
    data_frame_string = 'gene'

  r.assign('gene', expn_value)
  r('data = data.frame(' + data_frame_string + ')')
  coxuh_output = r('summary( coxph(formula = Surv(time, censor) ~ ' + formula_string + ', ' +
    'data = data, model=FALSE, x=FALSE, y=FALSE))')

  coef_ind = list(coxuh_output.names).index('coefficients')
  coeffs = coxuh_output[coef_ind]

  cox_dict = { 'name': gene_name }
  for multivariate in coeffs.rownames:
    cox_dict[multivariate] = {
      'z': coeffs.rx(multivariate, 'z')[0],
      'p': coeffs.rx(multivariate, 'Pr(>|z|)')[0]
    }

  return cox_dict

def write_file_with_results(input_file_name, results, outfile_location):
  input_file_name_slug = os.path.basename(input_file_name).split('.')[0]
  output_name = os.path.join(outfile_location, input_file_name_slug+'.out.csv')
  multivariates = results[0].keys()
  multivariates.remove('name')
  multivariates.remove('gene')
  multivariates.insert(0, 'gene')

  with open(output_name, 'w') as outfile:
    outfile.write('Gene/Probe,'+ ', '.join([m + ' Z Score, ' + m + ' P Value' for m in multivariates]) + '\n')
    for result in results:
      outfile.write(result['name'])
      for m in multivariates:
        outfile.write(
          ', ' +
          '{:g}'.format(result[m]['z']) + ', ' +
          '{:g}'.format(result[m]['p'])
        )
      outfile.write('\n')


def do_files(files, outdir, multivariates=[]):
  for f in files:
    try:
      survival_time, survival_censor, gene_names, patient_values, feature_names, features =  import_file(f)

      # select the requested features
      features_for_calculation = []
      feature_names_for_calculation = []
      for variable in multivariates:
        safe_feature_name = safe_string(variable)
        if safe_feature_name == 'all':
          if len(multivariates) == 1:
            features_for_calculation = features
          else:
            print 'Error: All features requested for multivariate calculation, but additional provided, ' + ', '.join(multivariates)
            sys.exit(2)
        if not safe_feature_name in feature_names:
          print 'Error: Requested multivariate ' + safe_feature_name + ' not found in given features, ' + ', '.join(features)
          sys.exit(2)
        else:
          feature_idx = feature_names.index(safe_feature_name)
          feature_names_for_calculation.append(safe_feature_name)
          features_for_calculation.append(features[feature_idx])

      results = []
      for i in range(len(patient_values)):
        results.append(coxuh(gene_names[i], patient_values[i] , survival_time , survival_censor, feature_names_for_calculation, features_for_calculation))
    except Exception as e:
      print "Something went wrong"
    finally:
      write_file_with_results(f, results, outdir)

def main(argv=None):
  if argv is None:
    argv = sys.argv
    try:
      opts, args = getopt.getopt(argv[1:], 'ho:i:vm:',
      ['help', 'input=', 'output=', 'multivariates='])
      ['help', 'input=', 'output=', 'multivariates=', 'interactive'])
    except getopt.error, msg:
      usage()

    infile = None
    outdir = '.'
    multivariates = []
    for option, value in opts:
      if option == '-v':
        verbose = True
      if option in ('-h', '--help'):
        usage()
      if option in ('-i', '--input'):
        infile = value
      if option in ('-o', '--output-directory'):
        outdir = value
      if option in ('-m', '--mutlivariates'):
        multivariates = value.split(',')

    if ('--interactive', '') in opts:
      try:
        survival_time, survival_censor, gene_names, patient_values, feature_names, features =  import_file_interactive(infile)
        print "File Imported Successfully, calculating survival..."
        results = []
        for i in range(len(patient_values)):
          results.append(coxuh(gene_names[i], patient_values[i] , survival_time , survival_censor, feature_names, features))
      except Exception as e:
        print "Something went wrong, ", e
      finally:
        print "Saving file..."
        write_file_with_results(infile, results, outdir)
        exit(0)

    if not infile:
      usage()
    do_files([infile], outdir, multivariates)

if __name__ == '__main__':
  main()

