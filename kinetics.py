from __future__ import division

import sys
import string
import math

from utils import ordered_dict, PrintObject
import nupack_out_grammar as ngram
from DNAfold import DNAfold
from multistrand import DNAkinfold

def read_nupack(filename):
  """Extracts the designed sequences and the mfe structures"""
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

class Complex(PrintObject):
  def __init__(self, strands, struct):
    self.strands = strands; self.struct = struct

def convert_ins(compl):
  """Convert internal representation for inputs into the simpler notation passed to kinfold."""
  strands = [strand.name for strand in compl.strands]
  return Complex(strands, compl.mfe_struct)   # TODO: should we use compl.struct instead?

def convert_outs(compl):
  """Convert output notation, much like convert_ins except that we run untill "DISASSOC"."""
  strands = [strand.name for strand in compl.strands]
  return Complex(strands, "DISASSOC")

def test_kinetics(kin, gate, trials=24, time=100000, temp=25, conc=10., num_proc=1, out_interval=-1):
  """Test times for inputs to combine/seperate into outputs"""
  ins  = [convert_ins(struct) for struct in kin.inputs]
  outs = [convert_outs(struct) for struct in kin.outputs]
  # HACK: We don't want to run until they've reached mfe structure.
  # Instead we wait for the strands to be in the right structures.
  # NOTE: This only tests that the 1st output structure is produced!
  outs_hack = outs[0:1]
  
  # Used strands is a set of all strands used in kinetic test, their names and sequences.
  used_strands = ordered_dict()
  for struct in ins:
    for strand_name in struct.strands:
      if strand_name not in used_strands:
        used_strands[strand_name] = gate.strands[strand_name].seq
  
  trials_each = int(math.ceil(trials / num_proc))
  return DNAkinfold(used_strands, ins, outs_hack, trials_each, time, temp, conc, num_proc, out_interval)
