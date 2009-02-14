"""Grammar for Shawn Ligocki's extended dot-paren notation"""

import sys

sys.path += ("..",) # Extend path to find pyparsing.py
from pyparsing import *

Map = lambda f: (lambda s,t,l: map(f, l))

int_ = Word(nums).setParseAction(Map(int))

symbol = Word(".()+", exact = 1)  # i.e. One of these symbols.
term = Optional(int_, default = 1) + symbol

expr = ZeroOrMore(Group(term))
# For Example: 6. 7( 4. + 7) 3. -> [[6, "."], [7, "("], [4, "."], [1, "+"], [7, ")"], [3, "."]]

def parse(s):
  """Parse a string in extended dot-paren notation."""
  return expr.parseString(s, parseAll=True)
