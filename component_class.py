import sys
import string

from utils import ordered_dict, PrintObject, error
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
  
  def assert_(self, statement, message):
    """Raise error if statement is false. Adds component information to message."""
    prefix = "In component %s: " % self.name
    if not statement:
      print error(prefix + message)
      sys.exit(1)
  
  ## Add information from document statements to object
  def add_sequence(self, name, const, length):
    if DEBUG: print "sequence", name
    self.assert_( name not in self.seqs, "Duplicate sequence definition for '%s'" % name )
    try:
      seq = Sequence(name, length, *const)
    except AssertionError, e:
      self.assert_(False, str(e))
    seq.full_name = self.prefix + name
    self.base_seqs[name] = self.seqs[name] = seq
  
  def add_super_sequence(self, name, const, length):
    if DEBUG: print "sup-sequence", name
    self.assert_( name not in self.seqs, "Duplicate sequence definition for '%s'" % name )
    for n, item in enumerate(const):
      if item[0] == "Sequence":
        seq_name, wc = item[1]
        if wc:
          const[n] = ~self.seqs[seq_name]
        else:
          const[n] =  self.seqs[seq_name]
    try:
      seq = SuperSequence(name, length, *const)
    except AssertionError, e:
      self.assert_(False, str(e))
    seq.full_name = self.prefix + name
    self.sup_seqs[name] = self.seqs[name] = seq
    # Add anonymous sequences
    for sub_seq in self.sup_seqs[name].seqs:
      if isinstance(sub_seq, AnonymousSequence):
        self.seqs[seq.name] = sub_seq
        self.base_seqs[seq.name] = sub_seq
  
  def add_strand(self, dummy, name, const, length):
    if DEBUG: print "strand", name
    self.assert_( name not in self.strands, "Duplicate strand definition for '%s'" % name )
    for n, item in enumerate(const):
      if item[0] == "Sequence":
        seq_name, wc = item[1]
        self.assert_( seq_name in self.seqs, "Sequence '%s' referenced before definion (in strand '%s')" % (seq_name, name) )
        if wc:
          const[n] = ~self.seqs[seq_name]
        else:
          const[n] =  self.seqs[seq_name]
    try:
      self.strands[name] = Strand(name, dummy, length, *const)
    except AssertionError, e:
      self.assert_(False, str(e))
    self.strands[name].full_name = self.prefix + name
    # Add anonymous sequences
    for seq in self.strands[name].seqs:
      if isinstance(seq, AnonymousSequence):
        self.seqs[seq.name] = seq
        self.base_seqs[seq.name] = seq
  
  def add_structure(self, opt, name, strands, struct):
    if DEBUG: print "struct", name
    self.assert_( name not in self.structs, "Duplicate structure definition for '%s'" % name )
    
    # Convert from list of strand names to list of strands
    for n, strand in enumerate(strands):
      self.assert_( strand in self.strands, "Strand '%s' referenced before definion (in structure '%s')" % (strand, name) )
      strands[n] = self.strands[strand]
    
    isdomain, struct = struct
    if isdomain: # This is a domain-based structure
      sub_structs = struct.split("+")
      full_struct = ""
      self.assert_( len(sub_structs) == len(strands), "Mismatched number of strands: Structure %s has %d strands, but structures %s implies %d." % (name, len(strands), struct, len(sub_structs)) )
      # For each strand expand out the structure
      for sub_struct, strand in zip(sub_structs, strands):
        self.assert_( len(sub_struct) == len(strand.seqs), "Mismatch: Strand %s in structure %s has %d domain(s), but sub-structure %s implies %d" % (strand.name, name, len(strand.seqs), sub_struct, len(sub_struct)) )
        for dp, domain in zip(sub_struct, strand.seqs):
          full_struct += dp * domain.length
        full_struct += "+"
      struct = full_struct[:-1] # Get rid of trailing +
    try:
      self.structs[name] = Structure(name, opt, struct, *strands)
    except AssertionError, e:
      self.assert_(False, str(e))
    self.structs[name].full_name = self.prefix + name
  
  def add_kinetics(self, inputs, outputs):
    if DEBUG: print "kin", self.kin_num
    for n, struct in enumerate(inputs):
      self.assert_( struct in self.structs, "Kinetic statement uses structure '%s' before it is defined." % struct )
      inputs[n] = self.structs[struct]
    for n, struct in enumerate(outputs):
      self.assert_( struct in self.structs, "Kinetic statement uses structure '%s' before it is defined." % struct )
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
      self.assert_( seq_name in self.seqs, "Declare statement references undefined sequence '%s'" % seq_name )
      if wc:
        self.input_seqs.append( self.seqs[seq_name].wc )
      else:
        self.input_seqs.append( self.seqs[seq_name] )
      
      if struct_name:
        self.assert_( struct_name in self.structs, "Declare statement references undefined structure '%s'" %  struct_name )
        self.input_structs.append( self.structs[struct_name] )
      else:
        self.input_structs.append(None)
    
    self.output_seqs = []
    self.output_structs = []
    for (seq_name, wc), struct_name in outputs:
      self.assert_( seq_name in self.seqs, "Declare statement references undefined sequence '%s'" % seq_name )
      if wc:
        self.output_seqs.append( self.seqs[seq_name].wc )
      else:
        self.output_seqs.append( self.seqs[seq_name] )
      
      if struct_name:
        self.assert_( struct_name in self.structs, "Declare statement references undefined structure '%s'" %  struct_name )
        self.output_structs.append( self.structs[struct_name] )
      else:
        self.output_structs.append(None)
        
  
  
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
