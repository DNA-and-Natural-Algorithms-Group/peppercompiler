import string

from utils import ordered_dict, PrintObject
from DNA_classes import *

DEBUG = False

class Gate(PrintObject):
  def __init__(self, name, params, inputs, outputs):
    """Initialized the gate with the declare statement"""
    self.decl_name = name
    self.params = params
    self.inputs  = [tuple(x) for x in inputs]
    self.outputs = [tuple(x) for x in outputs]
    
    self.seqs = ordered_dict()
    self.base_seqs = ordered_dict()  # Not super-sequences
    self.sup_seqs = ordered_dict()
    self.strands = ordered_dict()
    self.structs = ordered_dict()
    self.kinetics = ordered_dict()
    self.kin_num = 0
  
  ## Add information from document statements to object
  def add_sequence(self, name, const, length):
    if DEBUG: print "sequence", name
    assert name not in self.seqs, "Duplicate sequence definition"
    self.seqs[name] = Sequence(name, length, *const)
    self.base_seqs[name] = self.seqs[name]
  
  def add_super_sequence(self, name, const, length):
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
    # Add references junk sequences
    for seq in self.sup_seqs[name].seqs:
      if isinstance(seq, JunkSequence):
        self.seqs[seq.name] = seq
        self.base_seqs[seq.name] = seq
  
  def add_strand(self, dummy, name, const, length):
    if DEBUG: print "strand", name
    assert name not in self.strands, "Duplicate strand definition"
    for n, item in enumerate(const):
      if item[0] == "Sequence":
        if item[1][1] == "":
          const[n] =  self.seqs[item[1][0]]
        else:
          const[n] = ~self.seqs[item[1][0]]
    self.strands[name] = Strand(name, dummy, length, *const)
    # Add references junk sequences
    for seq in self.strands[name].seqs:
      if isinstance(seq, JunkSequence):
        self.seqs[seq.name] = seq
        self.base_seqs[seq.name] = seq
  
  def add_structure(self, opt, name, strands, struct):
    if DEBUG: print "struct", name
    assert name not in self.structs, "Duplicate structure definition"
    
    for n, strand in enumerate(strands):
      strands[n] = self.strands[strand]
    
    isdomain, struct = struct
    if isdomain: # This is a domain-based structure
      sub_structs = struct.split("+")
      struct = ""
      # For each strand expand out the structure
      for sub_struct, strand in zip(sub_structs, strands):
        assert len(sub_struct) == len(strand.seqs), (self.decl_name, name, strand.name, sub_struct, strand.seqs)
        for dp, domain in zip(sub_struct, strand.seqs):
          struct += dp * domain.length
        struct += "+"
      struct = struct[:-1] # Get rid of trailing +
    self.structs[name] = Structure(name, opt, struct, *strands)
  
  def add_kinetics(self, inputs, outputs):
    if DEBUG: print "kin", self.kin_num
    for n, struct in enumerate(inputs):
      assert struct in self.structs, "Kinetic statement in component '%s' uses structure '%s' before it is defined." % (self.decl_name, struct)
      inputs[n] = self.structs[struct]
    for n, struct in enumerate(outputs):
      assert struct in self.structs, "Kinetic statement in component '%s' uses structure '%s' before it is defined." % (self.decl_name, struct)
      outputs[n] = self.structs[struct]
    
    name = "Kin%d" % self.kin_num
    self.kin_num += 1
    self.kinetics[name] = Kinetics(name, list(inputs), list(outputs))
  
  
  def output_synthesis(self, prefix, outfile):
    """Output synthesis of all data into a single file."""
    if prefix:
      outfile.write("#\n## Gate %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top Gate\n")
    
    # Define sequences
    for seq in self.base_seqs.values():
      name = prefix + seq.name
      outfile.write("sequence %s = %s : %d\n" % (name, seq.const, seq.length))
    
    # Define super-sequences
    for seq in self.sup_seqs.values():
      name = prefix + seq.name
      const = string.join([seq.name for seq in seq.seqs], " ")
      outfile.write("sup-sequence %s = %s : %d\n" % (name, const, seq.length))
    
    # Define strands
    for strand in self.strands.values():
      name = prefix + strand.name
      const = string.join([seq.name for seq in strand.seqs], " ")
      if strand.dummy:
        dummy = "[dummy]"
      else:
        dummy = ""
      outfile.write("strands %s %s = %s : %d\n" % (dummy, name, const, seq.length))
    
    # Define structures
    for struct in self.structs.values():
      name = prefix + struct.name
      strands = string.join([strand.name for strand in struct.strands], " + ")
      outfile.write("structure [%dnt] %s = %s : %s\n" % (struct.opt, name, strands, struct.struct))
    
    # Define kinetics
    for kin in self.kinetics.values():
      # name = prefix + kin.name
      inputs = string.join([struct.name for struct in kin.inputs], " + ")
      outputs = string.join([struct.name for struct in kin.outputs], " + ")
      outfile.write("kinetic %s -> %s\n" % (inputs, outputs))
  
  
  def output_nupack(self, prefix, outfile):
    """Compile data into NUPACK format and output it"""
    if prefix:
      outfile.write("#\n## Gate %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top Gate\n")
    
    # Define structures
    for struct in self.structs.values():
      name = prefix + struct.name
      outfile.write("structure %s = %s\n" % (name, struct.struct))
      # TODO-maybe: test that all sequences are used.
    
    # Define sequences
    for seq in self.base_seqs.values():
      name = prefix + seq.name
      outfile.write("sequence %s = %s\n" % (name, seq.const))
    
    # Apply sequences to structures and set objective function
    for struct in self.structs.values():
      name = prefix + struct.name
      seqs = string.join([prefix+seq.name for seq in struct.base_seqs])
      outfile.write("%s : %s\n" % (name, seqs))
      if struct.opt: # Optimization parameter
        outfile.write("%s < %f\n" % (name, struct.opt))
    
    ### TODO: do something for prevents
