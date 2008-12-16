"""
The k-sequence avoidance algorithm attempts to construct sequences which avoid 
having any subsequences of length k which are complimentary unless explicitly
forced to be so.

Uses the input format from Winfree's SpuriousC algorithm.
"""

import random

import DNA_classes

def avoid(k, st, eq, wc):
  """Avoid k-subsequence complimentarity with constraints that nt be in st,
     and eq and wc specify representatives for equality and complimentarity."""
  st = [DNA_classes.group[x] for x in st]  # Get the groups from the letter.
  
