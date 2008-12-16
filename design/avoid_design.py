#!/usr/bin/env python
"""
Designs sequences using a k-sequence avoiding algorithm.
Uses Joe Zadah's input and output formats for compatibility with compiler.
"""

from spurious_design import prepare

def design(in_name, basename):
  prepare(in_name, basename)
  # TODO: run sequence avoider and process results

if __name__ == "__main__":
  import sys
  import re
  
  in_name = sys.argv[1]
  basename = re.sub(r"\.des\Z", "", in_name) # Makes *.des => *
  design(in_name, basename)
