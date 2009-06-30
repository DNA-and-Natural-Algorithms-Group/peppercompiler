import string

from utils import ordered_dict, PrintObject
from DNA_classes import *

DEBUG = False

class Component(PrintObject):
  """Stores information for a DNA component"""
  def __init__(self, name, prefix, params):
    self.name = name
    self.prefix = prefix  # Prefix for all sequences, structs, etc. in this component.
    self.params = params
    
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
    seq = Sequence(name, length, *const)
    seq.full_name = self.prefix + name
    self.base_seqs[name] = self.seqs[name] = seq
  
  def add_super_sequence(self, name, const, length):
    if DEBUG: print "sup-sequence", name
    assert name not in self.seqs, "Duplicate sequence definition"
    for n, item in enumerate(const):
      if item[0] == "Sequence":
        seq_name, wc = item[1]
        if wc:
          const[n] = ~self.seqs[seq_name]
        else:
          const[n] =  self.seqs[seq_name]
    seq = SuperSequence(name, length, *const)
    seq.full_name = self.prefix + name
    self.sup_seqs[name] = self.seqs[name] = seq
    # Add references junk sequences
    for sub_seq in self.sup_seqs[name].seqs:
      if isinstance(sub_seq, JunkSequence):
        self.seqs[seq.name] = sub_seq
        self.base_seqs[seq.name] = sub_seq
  
  def add_strand(self, dummy, name, const, length):
    if DEBUG: print "strand", name
    assert name not in self.strands, "Duplicate strand definition"
    for n, item in enumerate(const):
      if item[0] == "Sequence":
        seq_name, wc = item[1]
        if wc:
          const[n] = ~self.seqs[seq_name]
        else:
          const[n] =  self.seqs[seq_name]
    self.strands[name] = Strand(name, dummy, length, *const)
    self.strands[name].full_name = self.prefix + name
    # Add references junk sequences
    for seq in self.strands[name].seqs:
      if isinstance(seq, JunkSequence):
        self.seqs[seq.name] = seq
        self.base_seqs[seq.name] = seq
  
  def add_structure(self, opt, name, strands, struct):
    if DEBUG: print "struct", name
    assert name not in self.structs, "Duplicate structure definition for '%s'" % name
    
    # Convert from list of strand names to list of strands
    for n, strand in enumerate(strands):
      strands[n] = self.strands[strand]
    # TODO: strands = [self.strands[strand_name] for strand_name in strands]
    
    isdomain, struct = struct
    if isdomain: # This is a domain-based structure
      sub_structs = struct.split("+")
      full_struct = ""
      assert len(sub_structs) == len(strands), "Mismatched number of strands: structure %s.%s has %d strands, but structures %s implies %d." % (self.name, name, len(strands), struct, len(sub_structs))
      # For each strand expand out the structure
      for sub_struct, strand in zip(sub_structs, strands):
        assert len(sub_struct) == len(strand.seqs), "Mismatch: strand %s in structure %s.%s has %d domains, but sub-structure %s implies %d" % (strand.name, self.name, name, len(strand.seqs), sub_struct, len(sub_struct))
        for dp, domain in zip(sub_struct, strand.seqs):
          full_struct += dp * domain.length
        full_struct += "+"
      struct = full_struct[:-1] # Get rid of trailing +
    self.structs[name] = Structure(name, opt, struct, *strands)
    self.structs[name].full_name = self.prefix + name
  
  def add_kinetics(self, inputs, outputs):
    if DEBUG: print "kin", self.kin_num
    for n, struct in enumerate(inputs):
      assert struct in self.structs, "Kinetic statement in component '%s' uses structure '%s' before it is defined." % (self.name, struct)
      inputs[n] = self.structs[struct]
    for n, struct in enumerate(outputs):
      assert struct in self.structs, "Kinetic statement in component '%s' uses structure '%s' before it is defined." % (self.name, struct)
      outputs[n] = self.structs[struct]
    
    name = "Kin%d" % self.kin_num
    self.kin_num += 1
    self.kinetics[name] = Kinetics(name, list(inputs), list(outputs))
    self.kinetics[name].full_name = self.prefix + name
  
  def add_IO(self, inputs, outputs):
    """Add I/O information once we've read the component."""
    self.input_seqs = []
    self.input_structs = []
    for (seq_name, wc), struct_name in inputs:
      assert seq_name in self.seqs, "Declare statement in component '%s' references undefined sequence '%s'" % (self.name, seq_name)
      if wc:
        self.input_seqs.append( self.seqs[seq_name].wc )
      else:
        self.input_seqs.append( self.seqs[seq_name] )
      
      if struct_name:
        assert struct_name in self.structs, "Declare statement in component '%s' references undefined structure '%s'" % (self.name, struct_name)
        self.input_structs.append( self.structs[struct_name] )
      else:
        self.input_structs.append(None)
    
    self.output_seqs = []
    self.output_structs = []
    for (seq_name, wc), struct_name in outputs:
      assert seq_name in self.seqs, "Declare statement in component '%s' references undefined sequence '%s'" % (self.name, seq_name)
      if wc:
        self.output_seqs.append( self.seqs[seq_name].wc )
      else:
        self.output_seqs.append( self.seqs[seq_name] )
      
      if struct_name:
        assert struct_name in self.structs, "Declare statement in component '%s' references undefined structure '%s'" % (self.name, struct_name)
        self.output_structs.append( self.structs[struct_name] )
      else:
        self.input_structs.append(None)
        
  
  
  ## Outputs
  def output_synthesis(self, prefix, outfile):
    """Output synthesis of all data into a single file."""
    if prefix:
      outfile.write("#\n## Component %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top Component\n")
    
    # Define sequences
    for seq in self.base_seqs.values():
      name = prefix + seq.name
      outfile.write("sequence %s = %s : %d\n" % (name, seq.const, seq.length))
    
    # Define super-sequences
    for sup_seq in self.sup_seqs.values():
      name = prefix + sup_seq.name
      const = string.join([prefix + seq.name for seq in sup_seq.seqs], " ")
      outfile.write("sup-sequence %s = %s : %d\n" % (name, const, sup_seq.length))
    
    # Define strands
    for strand in self.strands.values():
      name = prefix + strand.name
      const = string.join([prefix + seq.name for seq in strand.seqs], " ")
      if strand.dummy:
        dummy = "[dummy]"
      else:
        dummy = ""
      outfile.write("strand %s %s = %s : %d\n" % (dummy, name, const, strand.length))
    
    # Define structures
    for struct in self.structs.values():
      name = prefix + struct.name
      strands = string.join([prefix + strand.name for strand in struct.strands], " + ")
      outfile.write("structure [%dnt] %s = %s : %s\n" % (struct.opt, name, strands, struct.struct))
    
    # Define kinetics
    for kin in self.kinetics.values():
      # name = prefix + kin.name
      inputs = string.join([prefix + struct.name for struct in kin.inputs], " + ")
      outputs = string.join([prefix + struct.name for struct in kin.outputs], " + ")
      outfile.write("kinetic %s -> %s\n" % (inputs, outputs))
  
  
  def output_nupack(self, prefix, outfile):
    """Compile data into NUPACK format and output it"""
    if prefix:
      outfile.write("#\n## Component %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top Component\n")
    
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
