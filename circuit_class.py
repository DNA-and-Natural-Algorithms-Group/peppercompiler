"""
The Circuit class stores all of the information in a circuit file and loads
  the gate/circuit temlates.
"""

import string

import DNA_classes
from utils import ordered_dict, PrintObject

DEBUG = False

def load_file(basename, args, path="."):
  """Load the file basename.(sys|comp) with args from path."""
  # Imported in the function to avoid circular import error.
  import os
  from circuit_parser import load_circuit
  from gate_parser import load_gate
  
  basename = os.path.join(path, basename) # Add the correct directory name
  new_path = os.path.dirname(basename) # Local directory of basename
  
  sys_name = basename+".sys"   # Name if file is a system specification
  comp_name = basename+".comp" # Name if file is a component spec.
  
  if os.path.isfile(sys_name):
    assert not os.path.isfile(comp_name), "Ambiguous specification: Both '%s' and '%s' exist. Please remove the one that does not belong and rerun the compiler." % (sys_name, comp_name)
    return load_circuit(sys_name, args, new_path)
  
  elif os.path.isfile(comp_name):
    return load_gate(comp_name, args)
  
  else:
    raise Exception, "Neither '%s' nor '%s' exist" % (sys_name, comp_name)

class Circuit(PrintObject):
  """Stores all the information in a circuit's connectivity file"""
  def __init__(self, path, name, params, inputs, outputs):
    """Initialized the cicuit with the declare statement"""
    self.path = path  # Local path to load gates relative to.
    
    self.decl_name = name
    self.inputs  = [tuple(x) for x in inputs]
    self.outputs = [tuple(x) for x in outputs]
    
    self.template = ordered_dict()
    self.glob = ordered_dict()
    self.lengths = ordered_dict()
    self.gates = ordered_dict()
    
    # Pointers to subgate's objects
    self.seqs = ordered_dict()
    self.nupack_seqs = ordered_dict()  # Not super-sequences
    self.sup_seqs = ordered_dict()
    self.strands = ordered_dict()
    self.structs = ordered_dict()
    self.kinetics = ordered_dict()
  
  ## Add information from document statements to object
  def add_import(self, *imports):
    if DEBUG: print "import", imports
    for path, name in imports:
      if name == None:
        # filename is used as the internal name by default
        if "/" not in path:
          name = path
        else:
          name = path[path.rfind("/")+1:]  # Strip off lower directories
      if DEBUG: print "import", name, path
      assert name not in self.template, "Duplicate import %s" % name
      self.template[name] = path

  def add_gate(self, gate_name, templ_name, templ_args, inputs, outputs):
    if DEBUG: print "gate", gate_name, templ_name, templ_args, inputs, outputs
    # Setup gates
    assert templ_name in self.template, "Template referenced before import: " + templ_name
    assert gate_name not in self.gates, "Duplicate gate definition: " + gate_name
    self.gates[gate_name] = this_gate = load_file(self.template[templ_name], templ_args, self.path)
    assert len(inputs) == len(this_gate.inputs),   "Length mismatch. %s / %s: %r != %r" % (gate_name, templ_name, len( inputs), len(this_gate.inputs ))
    assert len(outputs) == len(this_gate.outputs), "Length mismatch. %s / %s: %r != %r" % (gate_name, templ_name, len(outputs), len(this_gate.outputs))
    # Constrain all gate inputs and outputs
    ### TODO: marry these 2 together in a more eligent way.
    if isinstance(this_gate, Circuit): # If it's actually a circuit
      for (glob_name, glob_wc), (loc_name, loc_wc) in zip(list(inputs)+list(outputs), this_gate.inputs+this_gate.outputs):
        wc = (glob_wc != loc_wc)  # Are these signals complementary?
        if glob_name not in self.glob:
          self.glob[glob_name] = [(loc_name, gate_name, wc)]
          self.lengths[glob_name] = this_gate.lengths[loc_name]
        else:
          self.glob[glob_name].append( (loc_name, gate_name, wc) )
          assert self.lengths[glob_name] == this_gate.lengths[loc_name]
    else: # Otherwise it's a gate, so we want to constrain sequences
      for (glob_name, glob_wc), (loc_name, loc_wc) in zip(list(inputs)+list(outputs), this_gate.inputs+this_gate.outputs):
        wc = (glob_wc != loc_wc)  # Are these signals complementary?
        loc_seq = this_gate.seqs[loc_name]
        if glob_name not in self.glob:
          self.glob[glob_name] = [(loc_seq, gate_name, wc)]
          self.lengths[glob_name] = loc_seq.length
        else:
          self.glob[glob_name].append( (loc_seq, gate_name, wc) )
          assert self.lengths[glob_name] == loc_seq.length, (glob_name, self.lengths[glob_name], loc_seq.length)
    
    # Point to all objects in the gate
    # For each type of object: seqs, strands, ...
    for type_ in "seqs", "nupack_seqs", "sup_seqs", "strands", "structs", "kinetics":
      gate_objs = this_gate.__dict__[type_] # this_gate.seqs, this_gate.nupack_seqs, ...
      circuit_objs = self.__dict__[type_]   # self.seqs, self.nupack_seqs, ...
      # Point to all of those items from here with a prefix added to the name
      for name, obj in gate_objs.items():
        circuit_objs[gate_name + "-" + name] = obj
  
  def output_nupack(self, prefix, outfile):
    """Compile data into NUPACK format and output it"""
    if prefix:
      outfile.write("#\n## Subcircuit %s\n" % prefix[:-1])
    else:
      outfile.write("#\n## Top Circuit\n")
    # For each gate write it's contents in Zadeh's format with prefix "gatename-".
    for gate_name, template in self.gates.items():
      template.output_nupack(prefix+gate_name+"-", outfile)
    
    # For each global sequence connecting gates constrain them to be be equal.
    # To force this constraint I make  them all complimentary to a single dummy strand
    if prefix:
      outfile.write("#\n## Circuit %s Connectors\n" % prefix[:-1])
    else: 
      outfile.write("#\n## Top Circuit Connectors\n")
    for glob in self.glob:
      length = self.lengths[glob]
      glob_name = prefix + glob
      wc_name = glob_name + "-_WC"  # A dummy variable wc compliment to glob
      
      outfile.write("#\n## Global %s\n" % glob_name)
      outfile.write("sequence %s = %s\n" % (glob_name, "N" * length))
      outfile.write("sequence %s = %s\n" % (wc_name, "N" * length))
      
      dummy_name = "%s-_Self" % glob_name
      outfile.write("structure %s = %s\n" % (dummy_name, "(" * length + "+" + ")" * length))
      outfile.write("%s : %s %s\n" % (dummy_name, wc_name, glob_name))
      
      # For each instance of the global sequence, build a structure to 
      #   constrain it to the global one
      for loc_seq, gate_name, wc in self.glob[glob]:
        if isinstance(loc_seq, DNA_classes.Sequence):
          sig_name = gate_name + "-" + loc_seq.name
          seqs = prefix + sig_name
        elif isinstance(loc_seq, DNA_classes.SuperSequence):
          # If it's a super-seq, list the subsequences
          sig_name = gate_name + "-" + loc_seq.name
          seqs = string.join([prefix+gate_name+"-"+seq.name for seq in loc_seq.nupack_seqs])
        else: # it's a circuit signal
          sig_name = gate_name + "-" + loc_seq
          seqs = prefix + sig_name
        
        dummy_name = glob_name + "-" + sig_name
        outfile.write("structure %s = %s\n" % (dummy_name, "(" * length + "+" + ")" * length))
        
        if wc: # If it's complementary to the global signal, then we can enforce that directly
          outfile.write("%s : %s %s\n" % (dummy_name, glob_name, seqs))
        else:  # If it's equal, then we must force it complementary to the complement 'wc_name'
          outfile.write("%s : %s %s\n" % (dummy_name, wc_name, seqs))
