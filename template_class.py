import string
from generic_classes import ordered_dict, PrintObject
from DNA_classes2 import *

DEBUG = False

class Gate(PrintObject):
  def __init__(self):
    self.def_func = False
    self.seqs = ordered_dict()
    self.reg_seqs = ordered_dict()
    self.sup_seqs = ordered_dict()
    self.strands = ordered_dict()
    self.structs = ordered_dict()
    self.kinetics = ordered_dict()
    self.kin_num = 0
  
  ## Add information from document statements to object
  def add_function(self, (name, inputs, outputs)):
    if DEBUG: print "function", name
    assert not self.def_func, "Multiple function declarations"
    self.def_func = True
    self.func_name = name
    self.inputs = list(inputs)
    self.outputs = list(outputs)
  def add_sequence(self, (name, const, length)):
    if DEBUG: print "sequence", name
    assert not self.seqs.has_key(name), "Duplicate sequence definition"
    self.seqs[name] = Sequence(name, length, *const)
    self.reg_seqs[name] = self.seqs[name]
  def add_super_sequence(self, (name, const, length)):
    if DEBUG: print "sup-sequence", name
    assert not self.seqs.has_key(name), "Duplicate sequence definition"
    self.seqs[name] = SuperSequence(name, length, *const)
    self.sup_seqs[name] = self.seqs[name]
  def add_strand(self, (name, const, length)):
    if DEBUG: print "strand", name
    assert not self.strands.has_key(name), "Duplicate strand definition"
    self.strands[name] = Strand(name, length, *const)
  def add_structure(self, (no_mfe, name, strands, struct)):
    if DEBUG: print "struct", name
    assert not self.structs.has_key(name), "Duplicate structure definition"
    self.structs[name] = Structure(name, not no_mfe, struct, *strands)
  def add_kinetics(self, (inputs, outputs)):
    if DEBUG: print "kin", self.kin_num
    self.kinetics[self.kin_num] = Kinetics(self.kin_num, list(inputs), list(outputs))
    self.kin_num += 1
  
  def output_nupack(self, prefix, outfile):
    """Compile data into NUPACK format and output it"""
    # Define structures
    used_seqs = set()
    for struct in self.structs.values():
      name = prefix + struct.name
      outfile.write("structure %s = %s\n" % (name, struct.struct))
      used_seqs.update([ x for x in struct.nupack_seqs if not x.reversed] + \
                       [~x for x in struct.nupack_seqs if x.reversed])
    # Define sequences
    #raw_input(repr(used_seqs))
    for seq in used_seqs:
      #print seq.name
      name = prefix + seq.name
      outfile.write("sequence %s = %s\n" % (name, seq.nupack_constr))
    # Apply sequences to structures and set objective function
    for struct in self.structs.values():
      name = prefix + struct.name
      seqs = string.join([prefix+seq.name for seq in struct.nupack_seqs])
      outfile.write("%s : %s\n" % (name, seqs))
      if struct.mfe:
        outfile.write("%s < 1.0\n" % name)
    ### TODO: do something for prevents

