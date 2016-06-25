import pandas as pd
import glob
import os

def get_gpl_gse_from_path(f):
  sp = f.split('/')
  gpl = sp[1]
  gse_path = sp[2].split('.')[0]
  gse = gse_path.split('_')[0]
  return gpl, gse

path = os.path.join('multivariate-out', 'GPL*', 'GSE*')
files = glob.glob(path)

data_dict = {}
for f in files:
  gpl, gse = get_gpl_gse_from_path(f)
  print gpl, gse
  data = pd.read_csv(f, index_col=0, header=2)
  data_dict[gpl + '+' + gse] = data[u' gene Z Score']

print data_dict
combined_data = pd.DataFrame(data_dict)
combined_data.to_csv('combined_multivariate.csv')


