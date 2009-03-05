"""DNA design container classes"""
import string

# Global DNA nt groups
group = {"A": "A", "T": "T", "U": "T", "C": "C", "G": "G",
         "W": "AT", "S": "CG", "N": "ACGT"} #... Others can be put later if needed ...
rev_group = dict([(v, k) for (k, v) in group.items()])  # A reverse lookup for group.
complement = {"A": "T", "T": "A", "C": "G", "G": "C",
              "N": "N", "S": "S", "W": "W"} #... Others can be put later if needed ...
def wc(seq):
  """Returns the WC complement of a nucleotide sequence."""
  return string.join([complement[nt] for nt in reversed(seq)], "")

WILDCARD = "?"

class Sequence(object):
  """Container for sequences"""
  def __init__(self, name, length, *constraints):
    self.name = name
    self.length = length
    self.seq = None # Stores the sequence once it has been defined.
    self.reversed = False
    # Build the dummy sequence for the W-C complement
    self.wc = ReverseSequence(self)
    
    # Check length and resolve wildcard
    const = list(constraints)
    lengths = [num for num, code in const]
    wilds = lengths.count(WILDCARD)
    assert wilds in (0,1) , "Too many wildcards in sequence"
    if wilds == 0: # no wildcards
      check_length = sum(lengths)
      assert check_length == length, "Sequence length mismatch. %s: %r != %r" % (name, check_length, length)
    else: # one wildcard
      check_length = sum([x for x in lengths if x != WILDCARD])
      delta = length - check_length
      assert delta >= 0, "Sequence length mismatch"
      i = lengths.index(WILDCARD)
      const[i] = (delta, const[i][1])
    
    self.const = ""
    for (num, base) in const:
      self.const += base * num  # We represent constriants in long-form

  def __invert__(self):
    """Returns the Watson-Crick complementary sequence."""
    return self.wc
  def __repr__(self):
    return "Sequence(%(name)r, %(length)r, *%(const)r)" % self.__dict__

class ReverseSequence(Sequence):
  """Complements of defined sequences"""
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.length = wc.length
    self.seq = None # Stores the sequence once it has been defined.
    self.reversed = True
    self.wc = wc
  def __repr__(self):
    return "~Sequence(%(name)r, %(length)r, *%(const)r)" % self.wc.__dict__

class JunkSequence(Sequence):
  """Sequences we don't need lables for"""
  junk_num = 0
  def __init__(self, length, *const):
    name = "_Junk"+repr(JunkSequence.junk_num)
    Sequence.__init__(self, name, length, *const)
    JunkSequence.junk_num += 1


class SuperSequence(object):
  """Logical grouping of sequences"""
  def __init__(self, name, length, *constraints):
    self.name = name
    self.seqs = []
    self.seq = None # Stores the sequence once it has been defined.
    self.length = 0
    self.nupack_seqs = []
    wildcard = None
    # Process constraints
    for item in constraints:
      if isinstance(item, SuperSequence):
        # Expand previous super-sequence
        self.seqs.append(item)
        self.nupack_seqs += item.nupack_seqs
        self.length += item.length
      elif isinstance(item, Sequence):
        # If item is a previously defined sequence, leave it alone, it's good
        self.seqs.append(item)
        self.nupack_seqs.append(item)
        self.length += item.length
      else:
        # Otherwise it's a junk constraint
        ### TODO-maybe: combine adjacent junk sections together into one
        num, code = item
        if num != WILDCARD:
          junk_seq = JunkSequence(num, item)
          self.seqs.append(junk_seq)
          self.nupack_seqs.append(junk_seq)
          self.length += junk_seq.length
        else:
          assert not wildcard, "Multiple wildcards"
          wildcard = (len(self.seqs), len(self.nupack_seqs), item) # Index and entry of wildcard
          
    if not wildcard:
      assert self.length == length, "Super Sequence length mismatch. %s: %r != %r" % (name, self.length, length)
    else:
      delta = length - self.length
      assert delta >= 0, "Super Sequence too short. %s: %r > %r" % (name, self.length, length)
      i, j, item = wildcard
      junk_seq = JunkSequence(delta, item)
      self.seqs.insert(i, junk_seq)
      self.nupack_seqs.insert(j, junk_seq)
      self.length += junk_seq.length
    self.wc = ReverseSuperSequence(self)
  def __invert__(self):
    """Returns the Watson-Crick complementary sequence."""
    return self.wc
  def __repr__(self):
    return "SuperSequence(%(name)r, %(length)r, *%(seqs)r)" % self.__dict__

class ReverseSuperSequence(SuperSequence):
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.length = wc.length
    self.seqs = [~seq for seq in wc.seqs[::-1]]
    self.seq = None # Stores the sequence once it has been defined.
    self.nupack_seqs = [~seq for seq in wc.nupack_seqs[::-1]]
    self.wc = wc

class Strand(SuperSequence):
  """Container for strands"""
  def __init__(self, name, dummy, length, *constraints):
    SuperSequence.__init__(self, name, length, *constraints)
    self.dummy = dummy
  def __repr__(self):
    return "Strand(%(name)r, %(length)r, *%(seqs)r)" % self.__dict__

class Structure(object):
  """Container for structures/complexes"""
  def __init__(self, name, opt, struct, *strands):
    self.name = name
    self.opt = opt
    self.struct = struct
    self.mfe_struct = None # Stores the actual mfe structure once it's known
    self.strands = list(strands)
    self.seq = None # Stores the sequence once it has been defined.
    self.nupack_seqs = []
    strand_lengths = [len(strand) for strand in self.struct.split("+")] # Check that lengths match up
    for strand, length in zip(strands, strand_lengths):
      assert isinstance(strand, Strand), "Structure must get strands"
      assert strand.length == length, "Length mismatch"
      self.nupack_seqs += strand.nupack_seqs
  def __repr__(self):
    return "Structure(%(name)r, %(struct)r, *%(strands)r)" % self.__dict__
  def des_seqs(self):
    """Return the sequences list in .des format"""
    seqs = ""
    for seq in self.nupack_seqs:
      seqs += "%s " % seq.name
    return seqs

class Kinetics(object):
  def __init__(self, name, inputs, outputs):
    self.name = name
    self.inputs = inputs
    self.outputs = outputs      
  def __repr__(self):
    return "Kinetics(%(name)r, %(inputs)r, %(outputs)r)" % self.__dict__
