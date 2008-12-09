#!/usr/bin/env python
"""Designs a sequence randomly."""
from __future__ import division

import sys, random, string, re

# Extend path to see compiler library and the HU2dotParen library
here = sys.path[0] # System path to this module.
sys.path += (here+"/..", here+"/../HU2dotParen")
from DNAfold import DNAfold
from HU2dotParen import HU2dotParen
from nupack_in_parser import load_design
from DNA_classes import group, compliment

def random_choice(group):
  assert len(group) > 0
  return random.choice(group)

class default_list(list):
  def __init__(self, default_func):
    self.default = default_func
  def __getitem__(self, index):
    if len(self) <= index:
      self.extend([self.default(a) for a in range(len(self), index+1)])
    return list.__getitem__(self, index)

class cool_set(set):
  def __init__(self, *args):
    self.base = None
    self.const = None
    set.__init__(self, *args)

def set_str(string, i, val):
  assert len(val) == 1
  return string[:i] + val + string[i+1:]

class Connect(object):
  def __init__(self):
    self.data = default_list(lambda a: default_list(lambda b: [cool_set(), cool_set([(a,b)])]))
  def add(self, compl, x, y):
    # sequence number, location on sequence, parity (answers: is not wc compliment?)
    seq_x, loc_x, par_x = compl.seq_loc(x)
    seq_y, loc_y, par_y = compl.seq_loc(y)
    ### TODO: if seq_x, loc_x, par_x == seq_y, loc_y, par_y: return  # Speedup
    data_x = self.data[seq_x][loc_x]
    data_y = self.data[seq_y][loc_y]
    # Equate all involved appropriately
    data_x[par_x].update(data_y[not par_y])  
    data_x[not par_x].update(data_y[par_y])
    for i,j in data_x[True]:
      self.data[i][j][True]  = data_x[True]
      self.data[i][j][False] = data_x[False]
    for i,j in data_x[False]:
      self.data[i][j][True]  = data_x[False]
      self.data[i][j][False] = data_x[True]

def design(infilename, outfile):
  global d # DEBUG
  d = load_design(infilename)
  
  # Create complementarity matrix
  connect = Connect()
  for compl in d.complexes.values():
    for x,y in compl.bonds:
      connect.add(compl, x, y)
  # Randomly color sequences
  for i, seq in enumerate(d.seqs.values()):
    for j, symb in enumerate(seq.seq):
      data = connect.data[i][j]
      if not data[True].base:
        grp = set(group[symb])
        for (x,y) in data[True]:
          grp.intersection_update(group[d.seqs.get_index(x).seq[y]])
        for (x,y) in data[False]:
          grp.intersection_update( [compliment[symb] for symb in group[d.seqs.get_index(x).seq[y]]] )
        # Set constraints and choose random base
        data[True].const = tuple(grp)
        data[False].const = [compliment[symb] for symb in grp]
        data[True].base = random_choice(data[True].const)
        data[False].base = compliment[ data[True].base ]
      
      seq.seq = set_str(seq.seq, j, data[True].base)

  # Test thermo and redesign to improve structure
  build_structs(d.complexes)
  res = score(d.complexes)
  print len(res)
  good_enough = True
  while not good_enough:
    # Randomly mutate all bad structure
    for compl, index in res:
      seq_num, loc, par = compl.seq_loc(index)
      data = connect.data[seq_num][loc]
      data[True].base = random_choice(data[True].const)
      data[False].base = compliment[ data[True].base ]
    # Propogate mutations
    for i, seq in enumerate(d.seqs.values()):
      for j, symb in enumerate(seq.seq):
        seq.seq = set_str(seq.seq, j, connect.data[i][j][True].base)
    new_res = score(d.complexes)
    print len(new_res), len(res), len(best_res)
    res = new_res
    best_res = min(best_res, res)
    #raw_input()
    
  output_sequences(d, connect, outfile)

def build_structs(complexes):
  for compl in complexes.values():
    compl.dp_struct = HU2dotParen(compl.struct)
    breaks = [i for i, symb in enumerate(compl.dp_struct) if symb == "+"]
    lengths = [y - x - 1 for x, y in zip([-1] + breaks, breaks + [len(compl.dp_struct)])] # = diffs of breaks
    #print breaks, lengths
    lengths = iter(lengths)
    # Build strands list
    strand = []
    compl.strands = [strand]
    l = lengths.next()
    #print compl.seqs
    for seq in compl.seqs:
      if l == 0:
        l = lengths.next()
        strand = []
        compl.strands.append(strand)
      if seq.length <= l:
        strand.append(seq)
        #print compl.name, seq.name, l, type(l)
        l -= seq.length
      else:
        #print compl.name, seq.name, l, seq.length
        raise ValueError, "Sequences do not match structure"
    assert l == 0 and end_iter(lengths)

def end_iter(foo):
  try:
    foo.next()
    return False
  except StopIteration:
    return True

def build_seq(compl):
  return string.join([ string.join([seq.get_seq() for seq in strand], "") for strand in compl.strands], "+")

def score(complexes):
  """Find unwanted base-pairing"""
  diffs = []
  for compl in complexes.values():
    compl.seq = build_seq(compl)
    #print compl.name, compl.seq
    compl.mfe_struct, dG = DNAfold(compl.seq, 25)
    diffs += diff(compl, compl.dp_struct, compl.mfe_struct)
  return diffs

def diff(lable, str_a, str_b):
  ### TODO-maybe: make correct by having:
  ###       ..(((...+...)))..  to  ..(((..)+(..)))..  be 4 errors not just 2
  #print lable.name, str_a, str_b
  assert len(str_a) == len(str_b)
  str_a = str_a.replace("+", "")
  str_b = str_b.replace("+", "")
  return [(lable, i) for i, (a, b) in enumerate(zip(str_a, str_b)) if a != b]

def output_sequences(d, connect, fn):
  f = file(fn, "w")
  for compl in d.complexes.values():
    #print compl
    f.write("%d:%s\n" % (0, compl.name))
    gc_content = (compl.seq.count("C") + compl.seq.count("G")) / (len(compl.seq) - len(compl.strands))
    f.write("%s %f %f %d\n" % (compl.seq, 0, gc_content, 0))
    f.write("%s\n" % compl.dp_struct)
    f.write("%s\n" % compl.mfe_struct)
  for seq in d.seqs.values():
    # Write sequence (with dummy content)
    f.write("%d:%s\n" % (0, seq.name))
    gc_content = (seq.seq.count("C") + seq.seq.count("G")) / seq.length
    f.write("%s %f %f %d\n" % (seq.seq, 0, gc_content, 0))
    f.write(("."*seq.length+"\n")*2)
    # Write wc of sequence (with dummy content)
    seq = seq.wc
    f.write("%d:%s\n" % (0, seq.name))
    #wc_seq = string.join([compliment[symb] for symb in seq.seq[::-1]], "")
    f.write("%s %f %f %d\n" % (seq.get_seq(), 0, gc_content, 0))
    f.write(("."*seq.length+"\n")*2)
  f.write("Total n(s*) = %f" % 0)
  f.close()

if __name__ == "__main__":
  in_name = sys.argv[1]
  out_name = re.sub(r"\.des\Z", "", in_name) + ".mfe"
  design(in_name, out_name)

