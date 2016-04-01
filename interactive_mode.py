import sys
import help_message
import readline
import os

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

def repeating_row_selection(message, row_titles):
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

def repeating_gene_signature_row_set(prompt):
  selections = []
  readline.parse_and_bind("tab: complete")
  while True:
    probe_set_file = raw_input(prompt)
    if probe_set_file == 'END':
      break

    if not os.path.isfile(probe_set_file):
      print 'Error: File not found: ', probe_set_file
      sys.exit(1)
    with open(probe_set_file,'r') as f:
      probe_set = f.read().strip().split('\n')
    probe_set_name = os.path.basename(probe_set_file).split('.')[0]
    print probe_set_name
    print probe_set
    selections.append((probe_set_name, probe_set))
  return selections

def get_row_titles(name):
  with open(name, 'rU') as f:
    row_titles = []
    for i,line in enumerate(f.readlines()):
      row_title = line.split(',')[0]
      if row_title[0] == '!':
        continue
      if row_title == 'patient' or row_title == 'ID_REF':
        break
      elif i > 50:
        print "'patient' nor 'ID_REF' found. Aborting."
        help_message.usage()
      else:
        row_titles.append(row_title)
  return row_titles

def import_file_interactive(name):
    row_titles = get_row_titles(name)

    time_row_number = row_selection('Enter number for row that has survival time: ', row_titles)
    censor_row_number = row_selection('Enter number for row that has censor: ',row_titles)

    additional_variables_rows = []
    additional_variables_rows = repeating_row_selection('Enter a number for an additional variable, "END" to finish selection: ', row_titles)

    gene_signature_probe_sets = repeating_gene_signature_row_set('Enter a file containing a list of probes to calculate a gene signature as a multivariate, or "END" to finish selection: ')


    return {
        'name': name,
        'time_row_number': time_row_number,
        'censor_row_number': censor_row_number,
        'additional_variables_rows': additional_variables_rows,
        'gene_signature_probe_sets': gene_signature_probe_sets
        }
