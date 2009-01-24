from __future__ import division

import sys
import string

from utils import ordered_dict
import nupack_out_grammar as ngram
from DNAfold import DNAfold
from multi_kinfold import DNAkinfold

def read_nupack(filename):
  # Note: parseFile reads entire file into memory and then returns entire file of data
  stats, total_n_star = ngram.document.parseFile(filename)
  seqs = {}
  structs = {}
  # Catches everything (struct, seq and seq*) except bad inters (struct N struct)
  #print stats
  for stat in stats:
    #print stat
    name, seq, n_star, gc, mfe, ideal_struct, actual_struct = stat
    seqs[name] = seq
    structs[name] = actual_struct
  return seqs, structs

class Complex(object):
  def __init__(self, strands, struct):
    self.strands = strands; self.struct = struct

def test_kinetics(prefix, kin, seqs, mfe_structs, trials=24, time=100000, temp=25, conc=10., num_proc=1, out_interval=-1):
  """Test times for inputs to combine/seperate into outputs"""
  used_strands = ordered_dict()
  ## Subroutine
  def convert(foo):
    bar = []
    for compl in foo:
      # Load seq/struct with encoded name
      name = prefix + compl.name

      if seqs.has_key(name): # The standard easy method
        seq = seqs[name]
        struct = mfe_structs[name]
      else: # The fallback method
        seq = ""
        for strand in compl.strands:
          # For each strand combine sequences
          for small_seq in strand.nupack_seqs:
            seq += seqs[prefix + small_seq.name]
          seq += "+" # With strand break between strands
        seq = seq.rstrip("+")
        struct, dG = DNAfold(seq, temp) # We ignore the dG
        # Fill in sequence and mfe_structure info
        seqs[name] = seq
        mfe_structs[name] = struct

      # Load sequences for individual strands and check consistency
      strand_seqs = seq.split("+")
      assert len(compl.strands) == len(strand_seqs)
      these_strands = []
      for strand, strand_seq in zip(compl.strands, strand_seqs):
        these_strands.append(strand.name)
        assert len(strand_seq) == strand.length
        if not used_strands.has_key(strand.name):
          used_strands[strand.name] = strand_seq
        else:
          assert used_strands[strand.name] == strand_seq
        
      bar.append(Complex(these_strands, struct))
    return bar
  ## End Subroutine
  ins = convert(kin.inputs)
  outs = convert(kin.outputs)
  #HACK
  outs = outs[0:1]
  outs[0].struct = "DISASSOC"
  return DNAkinfold(used_strands, ins, outs, trials, time, temp, conc, num_proc, out_interval)


