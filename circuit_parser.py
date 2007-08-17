from circuit_class import Circuit
from pyparsing import *

## Some globals
K = CaselessKeyword
S = Suppress
O = Optional
List = lambda x: Group(ZeroOrMore(x))  # A grouped list
Map = lambda func: (lambda s,l,t: map(func, t) )  # A useful mapping function

import_ = "import"
gate = "gate"

ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
var_list = List(var + S(O("+")))  # Space sep list of variable names
path_chars = printables.replace(",", "")
path = Word(path_chars) # Path name in a directory structure
integer = Word(nums).setParseAction(Map(int))

# import Adder, HalfAdder5 as HalfAdder, templates/Crossing_Gates/LastAdder
import_stat = K(import_) + delimitedList(Group(path + O(S("as") + var, default=None)))
params = O( S("(") + Group(delimitedList(integer)) + S(")") , default=[])
gate_stat = K(gate) + var + S("=") + var + params + S(":") + var_list + S("->") + var_list

statement = import_stat | gate_stat

document = StringStart() + delimitedList(O(Group(statement)), delim="\n") + StringEnd()
document.ignore(pythonStyleComment)

def load_circuit(filename):
  """Load circuit connectivity file"""
  circuit = Circuit()
  
  # Build data
  statements = document.parseFile(filename)
  for stat in statements:
    #print list(stat)
    if stat[0] == import_:
      circuit.add_import(stat[1:])
    elif stat[0] == gate:
      circuit.add_gate(stat[1:])
    else:
      print stat
      raise Exception
  return circuit

