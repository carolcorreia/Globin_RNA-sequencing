#!/usr/bin/python3

__author__ = "Carolina Correia"
__email__ = "carolina.correia@ucdconnect.ie"
__date__ = "20th April 2017"

import sys
import os

# Declare local variables:
args = sys.argv
d = {}

# Check if all arguments were passed in the command line:
if len(args) != 2:
  print("You forgot to pass the file with IDs and new names.")
  sys.exit()

def get_dict():
  ''' 
  Read a two-column tab-delimited file and save the lines into
  a dictionary.

  '''
  with open(args[1], 'r') as file:
    for line in file:
      key, value = line.split()
      d[key] = value

  return d


get_dict()
for k in d:
  os.rename(k, d[k])