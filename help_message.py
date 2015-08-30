import sys

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
