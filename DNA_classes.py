"""DNA design container classes"""

WILDCARD = "?"

seqs = []
sup_seqs = []
strands = []
structs = []
kins = []

class Sequence(object):
  """Container for sequences"""
  def __init__(self, name, length, *constraints, **keys):
    self.name = name
    self.length = length
    self.const = list(constraints)
    self.reversed = False
    # Build the dummy sequence for the W-C complement
    self.wc = ReverseSequence(self)
    global seqs
    seqs.append(self)

    # Check length and resolve wildcard
    lengths = [num for num, code in constraints]
    wilds = lengths.count(WILDCARD)
    assert wilds in (0,1) , "Too many wildcards in sequence"
    if wilds == 0: # no wildcards
      check_length = sum(lengths)
      assert check_length == length, "Sequence length mismatch"
    else: # one wildcard
      check_length = sum([x for x in lengths if x != WILDCARD])
      delta = length - check_length
      assert delta >= 0, "Sequence length mismatch"
      i = lengths.index(WILDCARD)
      self.const[i] = (delta, self.const[i][1])

  def __invert__(self):
    """Returns the Watson-Crick complementary sequence."""
    return self.wc
  def __repr__(self):
    return "Sequence(%(name)r, %(length)r, %(const)r)" % self.__dict__
  def des_const(self):
    """Return constraints in .des format"""
    const = ""
    for num, code in self.const:
      const += "%d%s " % (num, code)
    return const

class ReverseSequence(Sequence):
  """Complements of defined sequences"""
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.length = wc.length
    self.reversed = True
    self.wc = wc
  def __repr__(self):
    return "~Sequence(%(name)r, %(length)r, %(const)r)" % self.wc.__dict__

class JunkSequence(Sequence):
  """Sequences we don't need lables for"""
  junk_num = 0
  def __init__(self, length, *const):
    #global junk_seqs
    #junk_seqs.append(self)
    name = "_Junk"+repr(JunkSequence.junk_num)
    Sequence.__init__(self, name, length, *const)
    JunkSequence.junk_num += 1


class SuperSequence(object):
  """Logical grouping of sequences"""
  def __init__(self, name, length, *constraints):
    self.name = name
    self.seqs = []
    self.length = 0
    wildcard = None
    global sup_seqs
    sup_seqs.append(self)
    # Process constraints
    for item in constraints:
      if isinstance(item, (Sequence, SuperSequence)):
        # If item is a previously defined sequence, leave it alone, it's good
        self.seqs.append(item)
        self.length += item.length
      else:
        # Otherwise it's a junk constraint
        ### TODO-maybe: combine adjacent junk sections together into one
        num, code = item
        if num != WILDCARD:
          junk_seq = JunkSequence(num, item)
          self.seqs.append(junk_seq)
          self.length += junk_seq.length
        else:
          assert not wildcard, "Multiple wildcards"
          wildcard = (len(self.seqs), item) # Index and entry of wildcard
          
    if not wildcard:
      assert self.length == length, "Super Sequence length mismatch"
    else:
      delta = length - self.length
      assert delta >= 0
      i, item = wildcard
      junk_seq = JunkSequence(delta, item)
      self.seqs.insert(i, junk_seq)
      self.length += junk_seq.length
    self.wc = None # lazy evaluate
  def __invert__(self):
    """Returns the Watson-Crick complementary sequence."""
    if not self.wc:
      self.wc = ReverseSuperSequence(self)
    return self.wc
  def __repr__(self):
    return "SuperSequence(%(name)r, %(length)r, %(seqs)r)" % self.__dict__
class ReverseSuperSequence(SuperSequence):
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.length = wc.length
    self.seqs = []
    for seq in wc.seqs[::-1]:
      self.seqs.append(~seq)
    self.wc = wc

class Strand(object):
  """Container for strands"""
  def __init__(self, name, length, *constraints):
    self.name = name
    self.seqs = []
    self.length = 0
    wildcard = None
    global strands
    strands.append(self)
    # Process constraints
    for item in constraints:
      if isinstance(item, (Sequence, SuperSequence)):
        # If item is a previously defined sequence, leave it alone, it's good
        self.seqs.append(item)
        self.length += item.length
      else:
        # Otherwise it's a junk constraint
        ### TODO-maybe: combine adjacent junk sections together into one
        num, code = item
        if num != WILDCARD:
          junk_seq = JunkSequence(num, item)
          self.seqs.append(junk_seq)
          self.length += junk_seq.length
        else:
          assert not wildcard, "Multiple wildcards"
          wildcard = (len(self.seqs), item) # Index and entry of wildcard
          
    if not wildcard:
      assert self.length == length, "Strand length mismatch"
    else:
      delta = length - self.length
      assert delta >= 0, "Strand '%s' too long '%d'" % (self.name, self.length)
      i, item = wildcard
      junk_seq = JunkSequence(delta, item)
      self.seqs.insert(i, junk_seq)
      self.length += junk_seq.length
  def __repr__(self):
    return "Strand(%(name)r, %(length)r, %(seqs)r)" % self.__dict__

class Structure(object):
  """Container for structures/complexes"""
  def __init__(self, name, struct, *strands):
    self.name = name
    self.strands = strands
    self.struct = struct
    global structs
    structs.append(self)
    for strand in strands:
      assert isinstance(strand, Strand), "Structure must get strands"
      ### TODO: check that structure lengths match strand lengths
  def __repr__(self):
    return "Structure(%(name)r, %(struct)r, %(strands)r)" % self.__dict__

class Kinetics(object):
  def __init__(self, name, inputs, outputs):
    self.name = name
    self.inputs = inputs
    self.outputs = outputs      
    global kins
    kins.append(self)
    # Make sure that both sides of the kinetic formula have the same strands
    in_strands = []
    out_strands = []
    for struct in inputs:
      assert isinstance(struct, Structure), "Kinetic must get structures, got (%s) of type (%s)" % (struct, type(struct))
      in_strands += struct.strands
    for struct in outputs:
      assert isinstance(struct, Structure), "Kinetic must get structures, got (%s) of type (%s)" % (struct, type(struct))
      out_strands += struct.strands
    for strand in in_strands:
      out_strands.remove(strand) # Assure all inputs are outputs (conservation of strands)
    assert not out_strands, "Kinetic must not output strands not input"
      
  def __repr__(self):
    return "Kinetics(%(name)r, %(inputs)r, %(outputs)r)" % self.__dict__

