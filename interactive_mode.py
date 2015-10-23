import sys
import numpy as np

def print_row_selection(row_titles):
  for i, title in enumerate(row_titles):
    print str(i+1) + '. ' + title

def convert_single_row_selection(row_selection, row_titles):
  try:
    row_selection_int = int(row_selection) - 1
    if row_selection_int >= len(row_titles):
      print 'Error: Invalid Selection'
      return None
    return row_selection_int
  except(ValueError):
    print 'Error: Invalid Selection'
    return None

def row_selection(message, row_titles, exit_on_invalid = True):
  print_row_selection(row_titles)
  row_selection = raw_input(message)
  row_selection = convert_single_row_selection(row_selection, row_titles)

  if row_selection == None:
    exit(1)
  return row_selection

def repeating_row_num_selection(message, row_titles):
  selections = []
  while True:
    print_row_selection(row_titles)

    selection = raw_input(message)
    if selection == 'END':
      break
    else:
      row_selection = convert_single_row_selection(selection, row_titles)
      if row_selection == None:
        continue
      selections.append(row_selection)
  return selections

def str_row_selection(selection_str, row_titles):
  if row_titles.count(selection_str) == 1:
    return row_titles.index(selection_str)

def repeating_string_row_selection(message, row_title_fn):
  selections = []
  row_titles = None

  while True:
    selection = raw_input(message)
    if selection == 'END':
      break
    else:
      if row_titles == None:
        row_titles = row_title_fn()
      row_selection = str_row_selection(selection, row_titles)
      if row_selection == None:
        continue
      print row_selection
      selections.append(row_selection)
  return selections

def get_row_titles(name):
  with open(name) as f:
    row_titles = []
    for line in f.readlines():
      row_title = line.split(',')[0]
      if row_title == 'patient':
        break
      else:
        row_titles.append(row_title)
  return row_titles

def import_file_interactive(name):
    row_titles = get_row_titles(name)

    time_row_number = row_selection('Enter number for row that has survival time: ', row_titles)
    censor_row_number = row_selection('Enter number for row that has censor: ',row_titles)

    additional_variables_rows = repeating_row_num_selection('Enter a number for an additional variable, "END" to finish selection: ', row_titles)

    print 'Next you can make gene selections. Note, after the first selection the gene names will be loaded'
    get_gene_rows = lambda: list(np.genfromtxt(name, usecols=0, delimiter=',', dtype=None, filling_values=''))
    gene_rows = repeating_string_row_selection('Enter a gene, "END" to finish selection: ', get_gene_rows)
    print gene_rows

    return name, time_row_number, censor_row_number, additional_variables_rows
