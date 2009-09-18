# Extend path to see compiler library
import sys
here = sys.path[0] # System path to this module.
sys.path.append(here+"/..")

from utils import ordered_dict, PrintObject

class Spec(PrintObject):
  def __init__(self):
    self.base_seqs = ordered_dict()
    self.sup_seqs = ordered_dict()
    self.seqs = ordered_dict() # All sequences
    self.strands = ordered_dict()
    self.structs = ordered_dict()
    self.equals = []
  
  def add_seq(self, name, template):
    assert name not in self.seqs, "Duplicate sequence definition: %s" % name
    self.base_seqs[name] = template # TODO: Make object, e.g. Sequence(name, template)
    self.seqs[name] = self.base_seqs[name]
  
  def add_sup_seq(self, name, sub_seq_names):
    assert name not in self.seqs, "Duplicate sequence definition: %s" % name
    self.sup_seqs[name] = sub_seq_names # TODO: Make object
    self.seqs[name] = self.sup_seqs[name]
  
  def add_strand(self, name, seq_names, dummy):
    assert name not in self.strands, "Duplicate strand definition: %s" % name
    self.strands[name] = seq_names, dummy # TODO: Make object
  
  def add_struct(self, name, strand_names, struct, params):
    assert name not in self.structs, "Duplicate structure definition: %s" % name
    self.strands[name] = strand_names, struct, params # TODO: Make object
  
  def add_equal(self, seq_names):
    self.equals.append( seq_names ) # TODO: Get objects?

