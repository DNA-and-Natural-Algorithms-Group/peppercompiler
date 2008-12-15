#!/usr/bin/env python
"""
Designs sequences using Winfree's SpuriousDesign/spuriousC.c algorithm.
Uses Joe Zadah's input and output formats for compatibility with compiler.
"""

from new_loading import load_file

def min_(foo):
  """Return the min element (or None if there are no elements)."""
  if not foo: # If foo is empty
    return None
  else:
    return min(foo)

def min_rep(a, b):
  """Return the minimum representative. Where None is no representative."""
  if a == None:
    return b
  elif b == None:
    return a
  else:
    return min(a, b)

class Connections(object):
  """
  Keeps track of which parts of dna sequence are constrained to be equal or 
  complimentary to which other parts.
  """
  def __init__(self, spec):
    # seq_const[n][0] = list of (s, i) where spec.structs[s].seq[m] is spec.seqs[n]
    # seq_const[n][1] = list of (s, i) where it's the compliment
    seq_const = [([], []) for seq in spec.seqs]
    
    self.struct_length = []
    # Build these sets.
    for s, struct in enumerate(spec.structs.values()):
      i = 0
      for seq in struct.seqs:
        # Add this struct, index to the appropriate set.
        seq_const[seq.num][seq.reversed].append((s, i))
        i += seq.length
      self.struct_length[s] = i
    
    self.table = {}
    # Expand these sets for specific nucleotides ...
    for seq_num, (eq, wc) in enumerate(seq_const):
      seq = spec.seqs.get_index(seq_num)
      for nt in xrange(seq.length):
        new_eq = [(s, i + nt) for (s, i) in eq])
        rep_eq = min_(new_eq)
        new_wc = [(s, i + (seq.length - 1 - nt)) for (s, i) in wc]
        rep_wc = min_(new_wc)
        # ... and save those into the connection table.
        for (s, i) in new_eq:
          self.table[(s, i)] = (rep_eq, rep_wc)
        for (s, i) in new_wc:
          self.table[(s, i)] = (rep_wc, rep_eq)
  
  def apply_wc(self, x, y):
    """Apply WC complimentarity condition between items."""
    eq_x, wc_x = self.table[x]
    eq_y, wc_y = self.table[y]
    # x is complimentary to y so merge their representatives.
    eq_x = wc_y = min_rep(eq_x, wc_y)
    eq_y = wc_x = min_rep(eq_y, wc_x)
    
    self.table[x] = eq_x, wc_x
    self.table[y] = eq_y, wc_y
  
  def printf(self):
    for z, (eq_z, wc_z) in self.table.items():
      print z
      print eq_z
      print wc_z
      print

def design(in_name, out_name):
  # Load specification from input file.
  spec = load_file(in_name)
  
  # Load structures and subsequences into Connections object
  # and connect sequences that are equal and complimentary.
  c = Connections(spec)
  #c.printf()
  # Connect regions that are constrained to helixes.
  for x, struct in enumerate(spec.structs.values()):
    for start, stop in struct.bonds:
      c.apply_wc((x, start), (x, stop))
  #c.printf()
  # Update the sequence constraints.
  
  # Convert from (strand, index) format to multi_index format
  
  # Create the conversion function
  f = {}
  j = 0
  for x, struct in enumerate(spec.structs.values()):
    for i in xrange(c.struct_length[s]):
      f[(x, i)] = j
      j += 1
    j += 2 # For double whitespace between structures
  
  table
      
  # Print SpuriousC style files
  

if __name__ == "__main__":
  import sys
  import re
  
  in_name = sys.argv[1]
  out_name = re.sub(r"\.des\Z", ".mfe", in_name) # Makes *.des => *.mfe
  design(in_name, out_name)


