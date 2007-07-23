"""DNA design container classes"""
import HU_parser

class Sequence(object):
  """Container for sequences"""
  def __init__(self, name, constraints, num):
    self.name = name
    self.constr = list(constraints)
    self.num = num
    self.reversed = False
    # Get sequence and length
    self.seq = ""
    for num, symb in constraints:
      self.seq += symb*num
    self.length = len(self.seq)
    # Build the dummy sequence for the W-C complement
    self.wc = ReverseSequence(self)
  def __invert__(self):
    """Returns the Watson-Crick complementary sequence."""
    return self.wc
  def __repr__(self):
    return "Sequence(%(name)r, %(constr)r)" % self.__dict__

class ReverseSequence(Sequence):
  """Complements of defined sequences"""
  def __init__(self, wc):
    self.name = wc.name + "*"
    self.length = wc.length
    self.num = wc.num
    self.reversed = True
    self.wc = wc
  def __repr__(self):
    return "~Sequence(%(name)r, %(constr)r)" % self.wc.__dict__

class Structure(object):
  """Container for structures/complexes"""
  def __init__(self, name, struct):
    self.name = name
    self.struct = struct
    self.bonds = HU_parser.get_bonds(struct)
  def set_seqs(self, seqs):
    assert not self.__dict__.has_key("seqs")
    self.seqs = tuple(seqs)

  def seq_loc(self, index):
    num = 0  # Current seq number
    while self.seqs[num].length  <=  index:
      index -= self.seqs[num].length  # Move past sequence
      num += 1

    if not self.seqs[num].reversed: # If it's normal direction
      return self.seqs[num].num, index, True
    else:  # It's backwards
      return self.seqs[num].num, self.seqs[num].length - index - 1, False

  def __repr__(self):
    return "Structure(%(name)r, %(struct)r)" % self.__dict__

