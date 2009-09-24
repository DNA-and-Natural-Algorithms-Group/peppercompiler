import string
from collections import defaultdict

from constraints import propagate_constraints
from PIL_parser import load_spec
from PIL_DNA_classes import group, rev_group, complement

def min_(xs):
  """Returns min element (or None if there are no elements)"""
  if len(xs) == 0:
    return None
  else:
    return min(xs)

def intersect_groups(x1, x2):
  g1 = group[x1]; g2 = group[x2]
  inter = set(g1).intersection(set(g2))
  # TODO: better errors for overconstrained system
  assert len(inter) > 0, "System over-constrained. %s and %s cannot be equal." % (x1, x2)
  inter = list(inter)
  inter.sort()
  inter = string.join(inter, "")
  return rev_group[inter]

class Constraints(object):
  def __init__(self):
    self.eq = {}
    self.wc = {}
    self.st = {}
  
  def init(self, x, letter="N"):
    """Must initialize an index 'x' before constraining it."""
    assert x not in self.eq and x not in self.wc and x not in self.st, "Cannot initialize an index twice: %r" % x
    self.eq[x] = []
    self.wc[x] = []
    self.st[x] = letter
  
  def add_eq(self, x, y):
    """Add a single equality constraint between items x and y."""
    self.eq[x].append(y)
    self.eq[y].append(x)
  
  def add_eqs(self, x, y, num):
    """Add equality constraints for 'num' bases starting with x and y"""
    for i in range(num):
      self.add_eq(x + i, y + i)
  
  def add_wc(self, x, y):
    """Add a single complementarity constraint between items x and y."""
    self.wc[x].append(y)
    self.wc[y].append(x)
  
  def add_wcs(self, x, y, num):
    """Add complementarity constraints for 'num' bases starting with x and y"""
    for i in range(num):
      self.add_wc(x + i, y - i)

  def propagate(self):
    self.eq, self.wc = propagate_constraints(self.eq, self.wc)
  
  def propagate_templates(self):
    """
    Propagate the constraints on sequence specified by the templates of equal or 
    complementary nucleotides.
    
    Warning: Only run after propagating the eq and wc constraints.
    """
    assert set(self.st.keys()) == set(self.eq.keys()) == set(self.wc.keys())
    done = set()
    for x, st in self.st.items():
      if x not in done:
        # Constraint must match all equal ...
        for y in self.eq[x]:
          assert y not in done, (x, y)
          st_y = self.st[y]
          st = intersect_groups(st, st_y)
        # ... and the complement of all wc
        for y in self.wc[x]:
          assert y not in done, (x, y)
          st_y = complement[self.st[y]]
          st = intersect_groups(st, st_y)
        
        # Apply new template
        for y in self.eq[x]:
          self.st[y] = st
          done.add(y)
        for y in self.wc[x]:
          self.st[y] = complement[st]
          done.add(y)
  
  def get_reps(self, isvalid=(lambda y: isinstance(y, int))):
    """Get representatives for eq and wc classes. Reps are minimum valid elements"""
    assert set(self.eq.keys()) == set(self.wc.keys())
    self.eq_rep = {}
    self.wc_rep = {}
    
    for x in self.eq.keys():
      self.eq_rep[x] = min_([y for y in self.eq[x] if isvalid(y)])
      self.wc_rep[x] = min_([y for y in self.wc[x] if isvalid(y)])
      for y in self.eq[x]:
        self.eq_rep[y] = self.eq_rep[x]
        self.wc_rep[y] = self.wc_rep[x]
      for y in self.wc[x]:
        self.eq_rep[y] = self.wc_rep[x]
        self.wc_rep[y] = self.eq_rep[x]
  
  def dump(self):
    """
    Output st constraints and eq and wc representatives as standard lists.
    These are not quite in spuriousC form because they still
      use None and base 0 indexing instead of 
      using -1 and base 1 indexing
    """
    assert set(self.st.keys()) == set(self.eq.keys()) == set(self.wc.keys())
    N = max([key for key in self.eq_rep.keys() if isinstance(key, int)]) + 1 # Number of indexes to be used
    eq = [self.eq_rep.get(i) for i in range(N)]
    wc = [self.wc_rep.get(i) for i in range(N)]
    st = [self.st.get(i) for i in range(N)]
    
    return eq, wc, st


def index_func_struct(spec):
  """Return an index function based on the STRUCTURE-centric specification"""
  struct_start = [None] * len(spec.structs)
  strand_start = [None] * len(spec.strands)
  for num, strand in enumerate(spec.strands.values()):
    strand.num = num
  
  prev_length = 0  # Keep a running sum of the lengths of structures
  for num, struct in enumerate(spec.structs.values()):
    struct.num = num
    struct_start[num] = prev_length
    for strand in struct.strands:
      if strand_start[strand.num] == None:
        strand_start[strand.num] = prev_length
      # On the linear constraint array we will put one blanks between strands
      prev_length += strand.length + 1
    # and two blanks between structures, 
    prev_length += 1
  
  def get_index(struct, index):
    """Return the general index for a nucleotide from the given structure."""
    assert index < struct.length, "get_index called with an index that is too big. %d >= %d" % (index, struct.length)
    result = struct_start[struct.num]
    for strand in struct.strands:
      # If this index puts us past this strand
      if index >= strand.length:
        index -= strand.length
        result += strand.length + 1
      # Otherwise it leaves us midway in this strand
      else:
        return result + index
    # We should never finish this loop (only happen if struct.length > sum(strands lengths)
    assert False, "Structure %s has length %d, but strands total length %d" % (struct.name, struct.length, sum([strand.length for strand in struct.strands]))
  
  def get_index_strand(strand, index):
    """Return the general index for a nucleotide from the given strand."""
    assert index < strand.length, "get_index called with an index that is too big. %d >= %d" % (index, strand.length)
    return strand_start[strand.num] + index
  
  return get_index, get_index_strand

def index_func_strand(spec):
  """Return an index function based on the STRAND-centric specification"""
  ### TODO: Deal with dummy strands
  strand_start = [None] * len(spec.strands)
  prev_length = 0  # Keep a running sum of the lengths of strands
  for num, strand in enumerate(spec.strands.values()):
    strand.num = num
    strand_start[num] = prev_length
    # On the linear constraint array we will put two blanks between strands
    prev_length += strand.length + 2
  
  def get_index_strand(strand, index):
    """Return the general index for a nucleotide from the given strand."""
    assert index < strand.length, "get_index called with an index that is too big. %d >= %d" % (index, strand.length)
    return strand_start[strand.num] + index
  
  def get_index(struct, index):
    """Return the general index for a nucleotide from the given structure."""
    assert index < struct.length, "get_index called with an index that is too big. %d >= %d" % (index, struct.length)
    for strand in struct.strands:
      # If this index puts us past this strand
      if index >= strand.length:
        index -= strand.length
      # Otherwise it leaves us midway in this strand
      else:
        return get_index_strand(strand, index)
    # We should never finish this loop (only happen if struct.length > sum(strands lengths)
    assert False, "Structure %s has length %d, but strands total length %d" % (struct.name, struct.length, sum([strand.length for strand in struct.strands]))
  
  return get_index, get_index_strand

def get_constraints_from_file(filename, struct_orient=False):
  """
  Load a file (in PIL format) and distil out the equality, complementarity 
  and templating constraints. 
  
  struct_orient == True if we are listing out all structures in constraint files
  struct_orient == False if we are listing out strands (only once each)
  """
  ## Load the specification
  spec = load_spec(filename)
  
  ## Create constraint object
  constraints = Constraints()
  
  
  ## Initialize all struct/strand constraints lists
  if struct_orient:
    get_index, get_index_strand = index_func_struct(spec)
    
    for struct in spec.structs.values():
      for x in range(struct.length):
        x2 = get_index(struct, x)
        constraints.init(x2)
    
    # Constrain all instances of the same strand to be equal
    for struct in spec.structs.values():
      offset = 0
      for strand in struct.strands:
        for x in range(strand.length):
          x2 = get_index_strand(strand, x)
          y2 = get_index(struct, offset + x)
          constraints.add_eq(x2, y2)
        offset += strand.length
  
  else: # if strand oriented
    get_index, get_index_strand = index_func_strand(spec)
    
    for strand in spec.strands.values():
      for x in range(strand.length):
        x2 = get_index_strand(strand, x)
        constraints.init(x2)
  
  
  ## Add constraints
  # Structural constraints
  for struct in spec.structs.values():
    for x, y in struct.bonds:
      x2 = get_index(struct, x)
      y2 = get_index(struct, y)
      constraints.add_wc(x2, y2)
  
  
  ## Initialize all sequence constraints lists
  num = 0
  for seq in spec.base_seqs.values():
    seq.num = num
    num += 1
    for x, letter in enumerate(seq.template):
      x2 = (seq.num, x)
      constraints.init(x2, letter)
    # We list both sequences and reverse_sequences making constraints code simpler
    seq.wc.num = num
    num += 1
    for x in range(seq.wc.length):
      x2 = (seq.wc.num, x)
      constraints.init(x2)
      constraints.add_wc(x2, (seq.num, seq.length - x - 1) ) # Add wc constraint
  # and for super-sequences
  for seq in spec.sup_seqs.values():
    seq.num = num
    num += 1
    for x in range(seq.length):
      x2 = (seq.num, x)
      constraints.init(x2)
    # and wc
    seq.wc.num = num
    num += 1
    for x in range(seq.wc.length):
      x2 = (seq.wc.num, x)
      constraints.init(x2)
      constraints.add_wc(x2, (seq.num, seq.length - x) ) # Add wc constraint
  
  
  # Equality constraints
  for eqlist in spec.equals:
    first_seq = eqlist[0]
    for seq in eqlist[1:]:
      assert seq.length == first_seq.length
      for x in range(seq.length):
        constraints.add_eq( (first_seq.num, x), (seq.num, x) )
  
  # Super-sequence constraints
  for seq in spec.sup_seqs.values():
    offset = 0
    for sub_seq in seq.seqs:
      for x in range(sub_seq.length):
        constraints.add_eq( (seq.num, offset + x), (sub_seq.num, x) )
      offset += sub_seq.length
  
  # Strand constraints
  for strand in spec.strands.values():
    offset = 0
    for seq in strand.seqs:
      for x in range(seq.length):
        x2 = get_index_strand(strand, offset + x)
        constraints.add_eq(x2, (seq.num, x) )
      offset += seq.length
  
  ## Propagate constraints
  constraints.propagate()
  
  ## Propagate templating constraints
  constraints.propagate_templates()
  
  ## Get constraints
  constraints.get_reps()
  eq, wc, st = constraints.dump()
  
  return eq, wc, st

if __name__ == "__main__":
  import sys
  
  print get_constraints_from_file(sys.argv[1], (len(sys.argv) > 2))
