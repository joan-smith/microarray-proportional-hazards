import zscores
import pandas as pd
import glob
import os
import sys
import errno
import numpy as np
from multiprocessing import Process


def make_feature_names(col):
  feature_names = []
  factor_string = 'factor{%s}: %s'
  for i, name in col.iteritems():
    if i in ['Time', 'Censor', 'GPL', 'GSE']:
      continue
    if i == 'Sex':
      feature_names.append(factor_string % ('M', name))
    elif i == 'Age':
      feature_names.append(name)
    elif i == 'Tumor Size':
      feature_names.append(name)
    else:
      feature_names.append(factor_string % ('0', name))

  return feature_names

def make_sure_path_exists(path):
  try:
    os.makedirs(path)
  except OSError as exception:
    if exception.errno != errno.EEXIST:
      raise

def do_threaded_work(multivariates):
  for i, col in multivariates.iteritems():
    gpl = col.loc['GPL']
    gse = col.loc['GSE'].upper()
    col = col.dropna()
    path = os.path.join('multivariate-clinical', 'Microarray-datasets', gpl, '*' + gse + '*.csv')
    fnames = glob.glob(path)
    print fnames
    if len(fnames) != 1:
      print "Failed to find file for: ", gpl, gse
      sys.exit(1)
    fname = fnames[0]
    feature_names = make_feature_names(col)
    out_path = os.path.join('multivariate-out', gpl)
    make_sure_path_exists(out_path)
    if len(glob.glob(os.path.join(out_path, gse + '*.csv'))) >= 1:
      continue
    zscores.script_run(fname,
      col.loc['Time'],
      col.loc['Censor'],
      metadata_feature_rows=feature_names,
      outdir=out_path)


def groupbytest(x):
  return x%4
def do_work():
  multivariates = pd.read_excel('multivariate-clinical/multivariate.xlsx', sheetname='Good', index_col=0, header=None)

  multivariates_splits = multivariates.groupby(groupbytest, axis=1)
  threads = []
  for x,y in multivariates_splits:
    print x
    p = Process(target=do_threaded_work, args=(y,))
    p.start()
    threads.append(p)

  for t in threads:
    t.join()


def main():
  do_work()

if __name__ == '__main__':
  main()

