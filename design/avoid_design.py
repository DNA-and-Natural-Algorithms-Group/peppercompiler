#!/usr/bin/env python
"""
Designs sequences using a k-sequence avoiding algorithm.
Uses Joe Zadah's input and output formats for compatibility with compiler.
"""

import avoid
from spurious_design import prepare

def design(in_name, k):
  st, eq, wc = prepare(in_name)
  d = avoid.Design(st, eq, wc)
  print d.avoid(k)
  # TODO: something with these sequences.

if __name__ == "__main__":
  import sys
  
  try:
    in_name = sys.argv[1]
    k = int(sys.argv[2])
  except:
    print "Usage: python avoid_design.py infilename k"
    sys.exit(1)
  
  design(in_name, k)

