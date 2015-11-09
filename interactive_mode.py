import sys
import help_message
import numpy as np
import readline

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

def get_all_row_titles(fname):
  return list(np.genfromtxt(fname, usecols=0, delimiter=',', dtype=None, filling_values=''))

def probe_file_selection(message, fname):
  file_continue = raw_input(message)
  if file_continue == 'no':
    return []

  readline.parse_and_bind("tab: complete")
  file_path = raw_input('Enter the path to the file with probes, one per line: ')
  try:
    probes_file = open(file_path, 'rU')
    probes = probes_file.readlines()
    probes = [probe.strip() for probe in probes]

    probe_rows = []
    row_titles = get_all_row_titles(fname)
    for probe in probes:
      if probe in row_titles:
        print 'found probe: ', probe
        probe_rows.append(row_titles.index(probe))
      else:
        print 'Warning: Probe', probe, 'not found in file.'

    return probe_rows
  except Exception as e:
    print 'Error: failed to read file provided at ' + file_path
    print e
    exit(1)

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
    censor_row_number = row_selection('Enter number for row that has censor: ',  row_titles)

    additional_variables_rows = []
    additional_variables_rows = repeating_row_selection('Enter a number for an additional variable, "END" to finish selection: ', row_titles)

    pcna25_rows = probe_file_selection('Add PCNA25? (yes or no): ', name)

    return {'name': name,
            'time_row': time_row_number,
            'censor_row': censor_row_number,
            'additional_variables': additional_variables_rows,
            'pcna25_rows': pcna25_rows }
