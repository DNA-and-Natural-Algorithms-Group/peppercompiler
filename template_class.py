import string
from generic_classes import ordered_dict, PrintObject
from DNA_classes import *

DEBUG = False

class Gate(PrintObject):
  def __init__(self, (name, params, inputs, outputs)):
    """Initialized the gate with the declare statement"""
    self.decl_name = name
    self.params = params
    self.in_sigs  = list(inputs)
    self.out_sigs = list(outputs)
    self.inputs  = [x for (x,y) in self.in_sigs]
    self.outputs = [x for (x,y) in self.out_sigs]
    
    self.seqs = ordered_dict()
    self.reg_seqs = ordered_dict()
    self.sup_seqs = ordered_dict()
    self.strands = ordered_dict()
    self.structs = ordered_dict()
    self.kinetics = ordered_dict()
    self.kin_num = 0
  
  ## Add information from document statements to object
  def add_sequence(self, (name, const, length)):
    if DEBUG: print "sequence", name
    assert name not in self.seqs, "Duplicate sequence definition"
    self.seqs[name] = Sequence(name, length, *const)
    self.reg_seqs[name] = self.seqs[name]
  
  def add_super_sequence(self, (name, const, length)):
    if DEBUG: print "sup-sequence", name
    assert name not in self.seqs, "Duplicate sequence definition"
    for n, item in enumerate(const):
      if item[0] == "Sequence":
        if item[1][1] == "":
          const[n] =  self.seqs[item[1][0]]
        else:
          const[n] = ~self.seqs[item[1][0]]
    self.seqs[name] = SuperSequence(name, length, *const)
    self.sup_seqs[name] = self.seqs[name]
  
  def add_strand(self, (name, const, length)):
    if DEBUG: print "strand", name
    assert name not in self.strands, "Duplicate strand definition"
    for n, item in enumerate(const):
      if item[0] == "Sequence":
        if item[1][1] == "":
          const[n] =  self.seqs[item[1][0]]
        else:
          const[n] = ~self.seqs[item[1][0]]
    self.strands[name] = Strand(name, length, *const)
  
  def add_structure(self, (mfe, name, strands, struct)):
    if DEBUG: print "struct", name
    assert name not in self.structs, "Duplicate structure definition"
    for n, strand in enumerate(strands):
      strands[n] = self.strands[strand]
    self.structs[name] = Structure(name, mfe, struct, *strands)
  
  def add_kinetics(self, (inputs, outputs)):
    if DEBUG: print "kin", self.kin_num
    for n, struct in enumerate(inputs):
      inputs[n] = self.structs[struct]
    for n, struct in enumerate(outputs):
      outputs[n] = self.structs[struct]
    self.kinetics[self.kin_num] = Kinetics(self.kin_num, list(inputs), list(outputs))
    self.kin_num += 1
  
  
  def output_nupack(self, prefix, outfile):
    """Compile data into NUPACK format and output it"""
    if prefix:
      outfile.write("#\n## Gate %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top Gate\n")
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
        outfile.write("%s < %f\n" % (name, struct.mfe))
    ### TODO: do something for prevents

