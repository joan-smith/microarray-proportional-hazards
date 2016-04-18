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
import rpy2.robjects as robjects
rpy2.robjects.numpy2ri.activate()

import help_message
import interactive_mode

def safe_string(string):
  safe_string = "".join([c for c in string if re.match(r'\w', c)])
  return safe_string

def clear_blanks(rows):
  for i, row in enumerate(rows):
    for j in range(len(row)):
      if row[j] == '':
        row[j] = np.nan
  return rows

def get_feature_row_names(patient_row_idx, cohort):
  return [name for name in cohort[0:patient_row_idx,0]]

def get_features(features_rows, patient_row_idx, cohort):
  if features_rows == None:
    feature_names = [name for name in cohort[2:patient_row_idx,0]]
    features = cohort[2:patient_row_idx,1:]
  else:
    feature_names = [name for name in cohort[features_rows, 0]]
    features = cohort[features_rows,1:]

  # replaces blanks with nans.
  features = clear_blanks(features)
  return features, feature_names

def get_formatted_row(row, formatter=lambda x: np.float(x)):
  formatted = []
  for i in row[1:]:
    if len(i) >= 1:
      formatted.append(formatter(i))
    else:
      formatted.append(np.nan)
  return formatted

def get_time_and_censor(cohort, time_row, censor_row):
  survival_time = get_formatted_row(cohort[time_row])
  try:
    survival_censor = get_formatted_row(cohort[censor_row], lambda x: np.int(x))
  except ValueError as e:
    print 'Unsupported Input Error: Time and censor provided must be integers. Row: ' + str(censor_row) + '.'
    print '   ', cohort[censor_row]
    print 'To work around, please use interactive mode and select supported time, censor, and covariates.'
    exit(1)

  return survival_time, survival_censor

def import_file(name, time_row=0, censor_row=1, features_rows=None):
  cohort = np.genfromtxt(name, delimiter=',', dtype=None, filling_values='', comments="!")
  survival_time, survival_censor = get_time_and_censor(cohort, time_row, censor_row)
  row_headers = list(cohort[:,0])

  if 'patient' in row_headers:
    patient_row_idx = row_headers.index('patient')
  elif 'ID_REF' in row_headers:
    patient_row_idx = row_headers.index('ID_REF')
  else:
    print 'Error: \'patient\' or \'ID_REF\' header not found'
    help_message.usage()

  # Note: these return *all* the metadata rows if feature_rows is None,
  # or the selected ones if feature_rows is set.
  features, feature_names = get_features(features_rows, patient_row_idx, cohort)
  all_metadata_row_names = get_feature_row_names(patient_row_idx, cohort)

  gene_names = row_headers[patient_row_idx+1:]
  patient_value = cohort[patient_row_idx+1:, 1:]
  patient_value = clear_blanks(patient_value)
  patient_value = patient_value.astype(np.float)

  input_data = {
      'patient_values': patient_value,
      'survival_time': survival_time,
      'survival_censor': survival_censor,
      'time_row_num': time_row,
      'censor_row_num': censor_row,
      'metadata_row_names': all_metadata_row_names,
      'gene_names': gene_names,
      'all_feature_names': feature_names,
      'all_features': features,
  }
  return input_data

def parse_gene_signatures(gene_names, patient_values, gene_signature_probe_sets):
  gene_signature_names = []
  gene_signatures = []
  #  for each probe:
  #    normalize probe values to 0 by average
  #    add normalized values to new array for probes
  #  average new array by column, treating nans appropriately.
  #  add averaged array to gene signature, and names to names
  for signature_name, gene_set in gene_signature_probe_sets:
    selected_genes = []
    normed_selected_gene_patient_values = np.empty((len(gene_set), patient_values.shape[1]))
    for i,gene in enumerate(gene_set):
      if not gene in gene_names:
        print 'Warn: probe/gene', gene, ' from ', signature_name, ' not found in input file'
        continue
      gene_index = gene_names.index(gene)
      gene_patient_values = patient_values[gene_index]
      # in order to ingore nans in the average, use nanmean
      avg_gene_patient_value = np.nanmean(gene_patient_values, dtype=np.float64)
      normed_gene_patient_values = np.subtract(gene_patient_values, avg_gene_patient_value)
      normed_selected_gene_patient_values[i] = normed_gene_patient_values

    averaged_gene_signature = np.nanmean(normed_selected_gene_patient_values, axis=0, dtype=np.float64)
    gene_signature_names.append(signature_name)
    gene_signatures.append(averaged_gene_signature)

  return gene_signature_names, gene_signatures

# params:
#  gene_name (string): name of gene to calculate multivariate for
#  expn_value (np array, float): expression value for every patient for that gene
#  surv_time (np array, float): survival time for each patient
#  surv_censor (np array, 0 or 1): survival censor for each patient
#  feature names (string list): list of names of multivariate features
#  features (np array float, shape=[len(feature_names)][len(patients)]): multivariate data for cox.
def coxuh(gene_name, expn_value, surv_time, surv_censor, feature_names, features):
  rpy.set_default_mode(rpy.NO_CONVERSION)
  r_old.library('survival')

  # remove missing data
  skip_cols = []
  for i in range(len(expn_value)):
    if np.isnan(expn_value[i]):
      skip_cols.append(i)

  if len(skip_cols) > (len(expn_value)/2):
    print 'warning: not enough samples for row ' + gene_name + ', skipping.'
    return {}

  expn_value = np.delete(expn_value, skip_cols)
  surv_time = np.delete(surv_time, skip_cols)
  surv_censor = np.delete(surv_censor, skip_cols)
  if len(feature_names) >= 1:
    features = np.delete(features, skip_cols, 1)


  r.assign('time', surv_time)
  r.assign('censor', surv_censor)
  safe_feature_names = []
  for idx, feature_name in enumerate(feature_names):
    if 'factor(' in feature_name:
      match =  re.search('factor\((.*)\): (.*)', feature_name)
      reference = match.group(1)
      factor_feature_name = safe_string(match.group(2))
      features = features[idx].astype(str)
      r.assign(factor_feature_name, robjects.FactorVector(features))
      # Once we have a feature set up in R, we need to set the reference level:
      #  r(feature_name <- relevel(feature_name, reference_level))
      r(factor_feature_name + ' <- relevel('+ factor_feature_name +', "' + reference + '")')
      safe_feature_names.append(factor_feature_name)
    else:
      feature = features[idx].astype(np.float)
      safe_feature_names.append(safe_string(feature_name))
      r.assign(safe_string(feature_name), feature)
  formula_string = ''
  if len(safe_feature_names) >= 1:
    formula_string = 'gene + ' + ' + '.join(safe_feature_names)
    data_frame_string = 'gene, '+ ', '.join(safe_feature_names)
  else:
    formula_string = 'gene'
    data_frame_string = 'gene'

  r.assign('gene', expn_value)
  r('data = data.frame(' + data_frame_string + ')')
  coxuh_output = r('summary( coxph(formula = Surv(time, censor) ~ ' + formula_string + ', ' +
    'data = data, model=FALSE, x=FALSE, y=FALSE))')

  coef_ind = list(coxuh_output.names).index('coefficients')
  coeffs = coxuh_output[coef_ind]

  patient_count_ind = list(coxuh_output.names).index('n')
  patient_count = coxuh_output[patient_count_ind][0]

  cox_dict = {
      'name': gene_name,
      'n': patient_count
  }
  for multivariate in coeffs.rownames:
    cox_dict[multivariate] = {
      'z': coeffs.rx(multivariate, 'z')[0],
      'p': coeffs.rx(multivariate, 'Pr(>|z|)')[0]
    }

  return cox_dict

def write_file_with_results(input_file_name, requested_data, results, outfile_location):
  input_file_name_slug = os.path.basename(input_file_name).split('.')[0]
  output_name = os.path.join(outfile_location, input_file_name_slug+'.out.csv')

  print "Writing file..."

  # get the list of requested multivariates to populate
  # column titles, and make sure that they're ordered well
  # with the gene first.
  multivariates = None
  i = 0
  while multivariates == None and i < len(results):
    if len(results[i].keys()) > 0:
      multivariates = results[i].keys()
      multivariates.remove('name')
      multivariates.remove('gene')
      multivariates.remove('n')
      multivariates.insert(0, 'gene')
    i += 1
  if i == len(results):
    print 'Error: no results'
    exit(1)


  time_row = requested_data['time_row_num']
  censor_row = requested_data['censor_row_num']
  with open(output_name, 'w') as outfile:
    outfile.write('Survival Time Row, ' + requested_data['metadata_row_names'][time_row] + ', ' + str(time_row+1) + ', Note: row number excludes rows beginning with "!" from row count' + '\n')
    outfile.write('Censor Row, ' + requested_data['metadata_row_names'][censor_row] + ', ' + str(censor_row+1) + '\n')
    print 'Gene/Probe, Patient Count, ' + ', '.join([m + ' Z Score, ' + m + ' P Value' for m in multivariates]) + '\n'
    outfile.write('Gene/Probe, Patient Count, ' + ', '.join([m + ' Z Score, ' + m + ' P Value' for m in multivariates]) + '\n')
    for result in results:
      if 'name' in result:
        outfile.write(result['name'])
        outfile.write(', {:d}'.format(result['n']))
        for m in multivariates:
          outfile.write(
            ', {:g}'.format(result[m]['z']) +
            ', {:g}'.format(result[m]['p'])
          )
        outfile.write('\n')
  print "Complete!"

def do_one_file(input_file, input_data, outdir="."):
  results = []
  gene_names = input_data['gene_names']
  patient_values = input_data['patient_values']
  survival_time = input_data['survival_time']
  survival_censor = input_data['survival_censor']
  feature_names = input_data['feature_names']
  features = input_data['features']
  try:
    for i in range(len(patient_values)):
      results.append(coxuh(gene_names[i], patient_values[i] , survival_time , survival_censor, feature_names, features))
  except Exception as e:
    print "Something went wrong"
    print e
  finally:
    write_file_with_results(input_file, input_data, results, outdir)

def do_files(files, outdir, multivariates=[]):
  for f in files:
    input_data = import_file(f)

    input_data['features'] = []
    input_data['feature_names'] = []
    # select the requested features
    for variable in multivariates:
      if variable == 'all':
        # make sure that user didn't ask for all + specific features
        if len(multivariates) == 1:
          input_data['features'] = input_data['all_features']
        else:
          print 'Error: All features requested for multivariate calculation, but additional provided, ' + ', '.join(multivariates)
          sys.exit(2)
      elif not variable in feature_names:
        print 'Error: Requested multivariate ' + variable + ' not found in given features, ' + ', '.join(features)
        sys.exit(2)
      else:
        feature_idx = feature_names.index(variable)
        input_data['feature_names'].append(variable)
        input_data['features'].append(input_data['all_features'][feature_idx])

      do_one_file(f, input_data, outdir)


def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'ho:i:vm:',
    ['help', 'input=', 'output-directory=', 'multivariates=', 'interactive'])
  except getopt.error, msg:
    help_message.usage()

  infile = None
  outdir = '.'
  multivariates = []
  for option, value in opts:
    if option == '-v':
      verbose = True
    if option in ('-h', '--help'):
      help_message.usage()
    if option in ('-i', '--input'):
      infile = value
    if option in ('-o', '--output-directory'):
      outdir = value
    if option in ('-m', '--mutlivariates'):
      multivariates = value.split(',')
  if not infile:
    help_message.usage()
  interactive = ('--interactive', '') in opts
  return infile, outdir, multivariates, interactive

def main(argv=None):
  if argv is None:
    argv = sys.argv
    infile, outdir, multivariates, interactive = get_options(argv)

    if interactive:
      selections = interactive_mode.import_file_interactive(infile)
      input_data = import_file(
          selections['name'],
          selections['time_row_number'],
          selections['censor_row_number'],
          features_rows=selections['additional_variables_rows'])

      gene_signature_names, gene_signatures = parse_gene_signatures(input_data['gene_names'], input_data['patient_values'], selections['gene_signature_probe_sets'])

      # in the interactive case, the features returned by import and
      # keyed by all_features, are the selected features.
      # so move them, and if necessary, add the new parsed gene_signature rows
      if len(gene_signature_names) == 0:
        multivariates = input_data['all_features']
        multivariate_names = input_data['all_feature_names']
      else:
        multivariates = np.vstack([gene_signatures, input_data['all_features']])
        multivariate_names = gene_signature_names + input_data['all_feature_names']
      input_data['features'] = multivariates
      input_data['feature_names'] = multivariate_names
      do_one_file(infile, input_data, outdir)
      sys.exit(0);
    else:
      do_files([infile], outdir, multivariates)

if __name__ == '__main__':
  main()

