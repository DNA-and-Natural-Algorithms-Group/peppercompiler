from circuit_class import Circuit
from pyparsing import *

## Some globals
K = CaselessKeyword
S = Suppress
O = Optional
List = lambda x: Group(ZeroOrMore(x))
import_ = "import"
input_ = "in"
gate = "gate"

ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
var_list = List(var + S(Optional("+")))  # Space sep list of variable names
path = Word(printables) # Path name in a directory structure

import_stat = S(K(import_)) + var
#long_import_stat = S(K(import_)) + path + S("as") var
input_stat  = S(K(input_))  + var
gate_stat = S(K(gate)) + var + S("=") + var + S(":") + var_list + S("->") + var_list

statement = import_stat | input_stat | gate_stat

document = StringStart() + delimitedList(O(statement), "\n") + StringEnd()
document.ignore(pythonStyleComment)

def load_circuit(filename):
  """Load circuit connectivity file"""
  circuit = Circuit()
  
  # Set parse actions
  import_stat.setParseAction(lambda s,t,l: circuit.add_import(l))
  input_stat.setParseAction(lambda s,t,l: circuit.add_input(l))
  gate_stat.setParseAction(lambda s,t,l: circuit.add_gate(l))
  
  # Build data
  document.parseFile(filename)
  return circuit

