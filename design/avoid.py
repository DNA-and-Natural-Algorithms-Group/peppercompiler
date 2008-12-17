"""
The k-sequence avoidance algorithm attempts to construct sequences which avoid 
having any subsequences of length k which are complimentary unless explicitly
forced to be so.

Uses the input format from Winfree's SpuriousC algorithm.
"""

import random

from spurious_design import NOTHING
import DNA_classes

def last(n, foo):
  """Gets the last n items in foo."""
  return foo[len(foo)-n:]  # Note that foo[-0:] would fail.

def avoid(k, st, eq, wc):
  """
  Avoid k-subsequence repeats or complementarity 
  with constraints that nt be in st, and 
  eq and wc specify representatives for equality and complimentarity.
  """
  assert len(st) == len(eq) == len(wc)
  # TODO: don't ignore breaks " "
  eq = [x for x, s in zip(eq, st) if s != " "]
  wc = [x for x, s in zip(wc, st) if s != " "]
  st = [DNA_classes.group[s] for s in st if s != " "]  # Get the groups from the letter.
  
  def step(i, part_seq, bad_seqs):
    """Step the algorithm at position i avoiding bad_sequences."""
    
    #print i, part_seq, bad_seqs
    
    # If we've assigned all bases, we're done.
    if i >= len(st):
      return part_seq
    
    # If we've already assigned an equal or wc constraint we must respect it
    if eq[i] < i:
      order = [part_seq[eq[i]]] # There is only one choice
    elif wc[i] != NOTHING and wc[i] < i:
      order = [DNA_classes.complement[part_seq[wc[i]]]] # There is only one choice
    
    # Otherwise, randomly order the allowed nts
    else:
      order = list(st[i])
      random.shuffle(order)
    
    #print order
    
    for nt in order:
      #print part_seq + nt, bad_seqs
      
      new_end = last((k-1), part_seq) + nt  # End of new seq with nucleotide added
      #print new_end
      # If the last 'k' nts are eq to a past seq ...
      if new_end in bad_seqs:
        j = bad_seqs[new_end]
        # ... and they aren't suposed to be, fail
        if eq[i-k:i] != eq[j-k:j]:
          #print "bad", part_seq + nt
          continue
      
      comp_end = DNA_classes.seq_comp(new_end)  # Complement
      #print comp_end
      # If the last 'k' nts are wc to a past seq ...
      if comp_end in bad_seqs:
        j = bad_seqs[comp_end]
        # ... and they aren't suposed to be, fail
        if wc[i-k:i] != eq[j-k:j]:
          #print "bad", part_seq + nt
          continue
      
      # Otherwise, use nt and step deeper.
      new_bad = bad_seqs.copy() # Don't mutate bad_seqs
      new_bad[new_end] = i
      
      # Try using this nt
      res = step(i+1, part_seq + nt, new_bad)
      if res:
        return res
      # Otherwise continue trying
      
      #print "worse", part_seq + nt
    
    return None # All nucleotides fail
  
  return step(0, "", {})

def test(k, n):
  """Try to find a k-sequence avoiding assignment for a length n single strand."""
  return avoid(k, "N"*n, range(n), [-1]*n)
