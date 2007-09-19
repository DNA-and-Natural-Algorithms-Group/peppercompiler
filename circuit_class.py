import string, time
from DNA_classes2 import Sequence
from generic_classes import ordered_dict, PrintObject
from template_parser import load_template

DEBUG = False

class Circuit(PrintObject):
  """Stores all the information in a circuit's connectivity file"""
  def __init__(self):
    self.template = ordered_dict()
    self.glob = ordered_dict()
    self.lengths = ordered_dict()
    self.gates = ordered_dict()
  
  ## Add information from document statements to object
  def add_import(self, imports):
    for path, name in imports:
	    if name == None:
	      # filename is used as the internal name by default
	      if "/" not in path:
	        name = path
	      else:
		      name = path[path.rfind("/")+1:]  # Strip off lower directories
	    if DEBUG: print "import", name, path
	    assert not self.template.has_key(name)
	    self.template[name] = path

  def add_gate(self, (gate_name, templ_name, templ_params, inputs, outputs)):
    if DEBUG: print "gate", gate_name, templ_name, templ_params, inputs, outputs
    # Setup gates
    self.gates[gate_name] = this_gate = load_template(self.template[templ_name], templ_params)
    assert len(inputs) == len(this_gate.inputs), "Length mismatch. %s / %s: %r != %r" % (gate_name, templ_name, len(inputs), len(this_gate.inputs))
    assert len(outputs) == len(this_gate.outputs)
    # Constrain all gate inputs and outputs
    for glob_name, loc_name in zip(list(inputs)+list(outputs), this_gate.inputs+this_gate.outputs):
      loc_seq = this_gate.seqs[loc_name]
      if not self.glob.has_key(glob_name):
        self.glob[glob_name] = [(loc_seq, gate_name)]
        self.lengths[glob_name] = loc_seq.length
      else:
        self.glob[glob_name].append((loc_seq, gate_name))
        assert self.lengths[glob_name] == loc_seq.length
  
  def output_nupack(self, filename):
    """Compile data into NUPACK format and output it"""
    outfile = file(filename, "w")
    outfile.write("## Specification compiled at %s\n" % time.ctime())
    # For each gate write it's contents in Zadeh's format with prefix "gatename-".
    for gate_name, template in self.gates.items():
      outfile.write("#\n## Gate %s\n" % gate_name)
      template.output_nupack(gate_name+"-", outfile)
    # For each global sequence connecting gates constrain them to be be equal.
    # To force this constraint I make  them all complimentary to a single dummy strand
    outfile.write("#\n#\n## Gate Connectors\n")
    for glob_name in self.glob:
      length = self.lengths[glob_name]
      outfile.write("#\n## Global %s\n" % glob_name)
      outfile.write("sequence %s = %dN\n" % (glob_name, length))
      for loc_seq, gate_name in self.glob[glob_name]:
        dummy_name = "%s-%s" % (glob_name, gate_name)
        outfile.write("structure %s = H%d(+)\n" % (dummy_name, length))
        if isinstance(loc_seq, Sequence):
          seqs = gate_name+"-"+loc_seq.name
        else:
          seqs = string.join([gate_name+"-"+seq.name for seq in loc_seq.nupack_seqs])
        outfile.write("%s : %s %s\n" % (dummy_name, glob_name, seqs))

