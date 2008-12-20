"""
The k-sequence avoidance algorithm attempts to construct sequences which avoid 
having any subsequences of length k which are complimentary unless explicitly
forced to be so.

Allows complementarity on overlapping regions.
Thus k=3 allows ACGT even though ACG ~ CGT

Uses the input format from Winfree's SpuriousC algorithm.
"""

import sys
import random

from spurious_design import NOTHING
import DNA_classes

def last(n, foo):
  """Gets the last n items in foo."""
  if n == 0:
    return foo[len(foo):]
  else:
    return foo[-n:]  # Note that foo[-0:] would fail.

def get_group(s):
  """Convert a letter into a group. i.e. "S" -> "CG", "N" -> "ATCG", etc. """
  if s == " ":
    return " "
  else:
    return DNA_classes.group[s]

def avoid(k, st, eq, wc):
  """
  Avoid k-subsequence repeats or complementarity 
  with constraints that nt be in st, and 
  eq and wc specify representatives for equality and complimentarity.
  """
  assert len(st) == len(eq) == len(wc)
  n = len(st)
  sys.setrecursionlimit(max(10*n, 1000)) # Might fail on to.dna for >25,000 recursion
  st = map(get_group, st)  # Get the groups from the letter.
  
  
  def step(i, part_seq, bad_seqs):
    """Step the algorithm at position i avoiding bad_sequences."""
    # If we've assigned all bases, we're done.
    if i >= len(st):
      return part_seq
    
    # If it's a strand break, pop it on and keep going.
    if st[i] == " ":
      return step(i+1, part_seq + " ", bad_seqs)
    
    # If we've already assigned an equal or wc constraint we must respect it
    if eq[i] < i:
      order = [part_seq[eq[i]]] # There is only one choice
    elif wc[i] != NOTHING and wc[i] < i:
      order = [DNA_classes.complement[part_seq[wc[i]]]] # There is only one choice
    
    # Otherwise, randomly order the allowed nts
    else:
      order = list(st[i])
      random.shuffle(order)

    # Try each nucleotide in group until one works
    for nt in order:
      new_end = last((k-1), part_seq) + nt  # End of new seq with nucleotide added
      # If the lask k nts are in a single strand
      if " " not in new_end:
        # If the last 'k' nts are eq to a past seq ...
        if new_end in bad_seqs:
          j = bad_seqs[new_end]
          # ... and they aren't suposed to be, fail
          if eq[i+1-k:i+1] != eq[j+1-k:j+1]:
            continue
        
        comp_end = DNA_classes.seq_comp(new_end)  # Complement
        # If the last 'k' nts are wc to a past seq ...
        if comp_end in bad_seqs:
          j = bad_seqs[comp_end]
          # ... and they aren't suposed to be, fail.
          # Note: we have reversed the indexing to wc.
          # TODO-test: ignore overlapping regions (j > i-k) because they will not bond.
          #if j <= i-k and wc[i:i-k:-1] != eq[j+1-k:j+1]:
          if wc[i:i-k:-1] != eq[j+1-k:j+1]:
            continue
      
      # Otherwise, use nt and step deeper.
      new_bad = bad_seqs.copy() # Don't mutate bad_seqs
      new_bad[new_end] = i
      
      # Try using this nt
      res = step(i+1, part_seq + nt, new_bad)
      if res:
        return res
      # Otherwise continue trying
      
    return None # All nucleotides fail
  
  return step(0, "", {})

def testU(k, n):
  """Try to find a k-sequence avoiding assignment for a length n unpaired single strand."""
  return avoid(k, "N"*n, range(n), [-1]*n)

def testH(k, n):
  """Try to find a k-sequence avoiding assignment for a length n helix."""
  return avoid(k, "N"*n + " " + "N"*n, range(2*n+1), range(2*n, -1, -1))
