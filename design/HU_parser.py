"""Grammar for Joe Zadeh's Helix Unpaired format"""

from pyparsing import *

Map = lambda f: (lambda s,t,l: map(f, l))

int_ = Word(nums).setParseAction(Map(int))

# This is a recursive grammar and requires using the 'exression' before it is defined
expr = Forward()

strand_break = Group("+")
unpaired = Group("U" + int_)
helix = Group("H" + int_ + Suppress("(") + Group(expr) + Suppress(")"))

term = strand_break | unpaired | helix

# Now we can define expression which is recursively defined
expr << ZeroOrMore(term)

statement = StringStart() + expr + StringEnd()

parse = statement.parseString

# For Example: U6 H7(U4 +) U3 -> [["U", 6], ["H", 7, [["U", 4], ["+"]]], ["U", 3]]

class stack(list):
  def push(self, val):
    self.append(val)
  def pop(self):
    return list.pop(self, -1)

class Bonds(object):
  def __init__(self):
    self.index = 0
    self.bonds = []

def get_bonds(struct):
  bonds = Bonds()
  ## Recursive function for finding bonds
  def recurse(p):
    for foo in p:
      # case foo = ["+"]
      if foo[0] == "+":
        continue
      # case foo = ["U", num]
      elif foo[0] == "U":
        bonds.index += foo[1]  # pass over unpaired 
      # case foo = ["H", num, sub_expr]
      else: # foo[0] == "H":
        start_index = bonds.index
        num_bonds = foo[1]
        bonds.index += num_bonds
        # Pass the helix to look inside
        recurse(foo[2])
        # Close the helix
        bonds.index += num_bonds
        end_index = bonds.index
        for n in range(num_bonds):
          bonds.bonds.append( (start_index + n, end_index-1 - n) )
  ## End of recursive function
  par = statement.parseString(struct)
  recurse(par)
  return bonds

