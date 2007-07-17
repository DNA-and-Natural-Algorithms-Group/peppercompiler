import string
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
  def add_import(self, (name,)):
    if DEBUG: print "import", name
    ### TODO: allow different name than path (i.e. import GeorgHalfAdder0712 as HalfAdder)
    path = name
    assert not self.template.has_key(name)
    self.template[name] = load_template(path)
  def add_input(self, (name,)):
    pass
  #  if DEBUG: print "in", name
  #  # Setup inputs (not strictly necessary, may scrap ...)
  #  self.glob[name] = []
  #  self.lengths[name] = None
  def add_gate(self, (gate_name, templ_name, inputs, outputs)):
    if DEBUG: print "gate", gate_name
    # Setup gates
    this_templ = self.template[templ_name]
    self.gates[gate_name] = this_templ
    assert len(inputs) == len(this_templ.inputs)
    assert len(outputs) == len(this_templ.outputs)
    # Constrain all gate inputs and outputs
    for glob_name, loc_name in zip(list(inputs)+list(outputs), this_templ.inputs+this_templ.outputs):
      loc_seq = this_templ.seqs[loc_name]
      if not self.glob.has_key(glob_name):
        self.glob[glob_name] = [(loc_seq, gate_name)]
        self.lengths[glob_name] = loc_seq.length
      else:
        self.glob[glob_name].append((loc_seq, gate_name))
        assert self.lengths[glob_name] == loc_seq.length
  
  def output_nupack(self, filename):
    """Compile data into NUPACK format and output it"""
    outfile = file(filename, "w")
    # For each gate write it's contents in Zadeh's format with prefix "gatename-".
    for gate_name, template in self.gates.items():
      outfile.write("## Gate %s\n" % gate_name)
      template.output_nupack(gate_name+"-", outfile)
    # For each global sequence connecting gates constrain them to be be equal.
    # To force this constraint I make  them all complimentary to a single dummy strand
    outfile.write("## Gate Connectors\n")
    for glob_name in self.glob:
      length = self.lengths[glob_name]
      outfile.write("## Global %s\n" % glob_name)
      outfile.write("sequence %s = %dN\n" % (glob_name, length))
      for loc_seq, gate_name in self.glob[glob_name]:
        dummy_name = "%s-%s" % (glob_name, gate_name)
        outfile.write("structure %s = H%d(+)\n" % (dummy_name, length))
        if isinstance(loc_seq, Sequence):
          seqs = loc_seq.name
        else:
          seqs = string.join([gate_name+"-"+seq.name for seq in loc_seq.nupack_seqs])
        outfile.write("%s : %s %s\n" % (dummy_name, glob_name, seqs))

