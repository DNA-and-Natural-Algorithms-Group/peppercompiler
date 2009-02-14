"""Grammar for Joe Zadeh's Helix Unpaired format"""

import sys
sys.path += ("..",) # Extend path to find pyparsing.py
from pyparsing import *

Map = lambda f: (lambda s,t,l: map(f, l))

# Container classes
class StrandBreak(object):
  def __init__(self, s, t, l):
    assert l[0] == "+"
class Unpaired(object):
  def __init__(self, s, t, l):
    assert l[0] == "U"
    self.length = l[1]
class Helix(object):
  def __init__(self, s, t, l):
    # l = [ ["H", len], contents]
    assert l[0][0] == "H"
    self.length = l[0][1]
    self.contents = list(l[1])

int_ = Word(nums).setParseAction(Map(int))

# This is a recursive grammar and requires using the 'exression' before it is defined
expr = Forward()

strand_break = Literal("+").setParseAction(StrandBreak)
unpaired = ("U" + int_).setParseAction(Unpaired)

helix_term = Group("H" + int_)
helix = (helix_term + Suppress("(") + Group(expr) + Suppress(")")).setParseAction(Helix)

term = strand_break | unpaired | helix

# Now we can define expression which is recursively defined
expr << ZeroOrMore(term)

# For Example: U6 H7(U4 +) U3 -> [["U", 6], [["H", 7], [["U", 4], "+"]], ["U", 3]]
