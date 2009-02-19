import sys

from circuit_class import Circuit
from gate_parser import substitute
from var_substitute import process

from pyparsing import *

## Some globals
K = CaselessKeyword
S = Suppress
O = Optional
List = lambda x: Group(ZeroOrMore(x))  # A grouped list
Map = lambda func: (lambda s,l,t: map(func, t) )  # A useful mapping function

decl = "declare"
system = "system"
import_ = "import"
gate = "component"

# Don't ignore newlines!
ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
var_list = List(var + S(O("+")))  # Space sep list of variable names
path = Word(alphanums+".-_/~") # Path name in a directory structure
py_chars = printables.replace(",", "").replace(")", "")
python_object = Word(py_chars, py_chars+" ").setParseAction(Map(eval))

# declare system <cicuit name> = <func name>(<params>): <inputs> -> <outputs>
decl_params = O(S("(") + Group(delimitedList(var)) + S(")"), default=[])
decl_stat = K(decl) + S(system) + var + decl_params + S(":") + var_list + S("->") + var_list
# import Adder, HalfAdder5 as HalfAdder, templates/Crossing_Gates/LastAdder
import_stat = K(import_) + delimitedList(Group(path + O(S("as") + var, default=None)))
# gate <name> = <template name>(<params>): <inputs> -> <outputs>
gate_params = O( S("(") + Group(delimitedList(python_object)) + S(")") , default=[])
gate_stat = K(gate) + var + S("=") + var + gate_params + S(":") + var_list + S("->") + var_list

statement = import_stat | gate_stat

document = StringStart() + Group(decl_stat) + S("\n") + \
           Group(delimitedList(O(Group(statement)), delim="\n")) + StringEnd()
document.ignore(pythonStyleComment)

def load_circuit(filename, args):
  """Load circuit connectivity file"""
  try:
    # Open file and do parameter substitution
    doc = substitute(filename, args)
  except ParseBaseException, e:
    print
    print "Parsing error in circuit:", filename
    print e
    sys.exit(1)
    
  try:
    # Load data
    decl_val, statements = document.parseString(doc)
  except ParseBaseException, e:
    print
    print doc
    print "Parsing error in circuit:", filename
    print e
    sys.exit(1)
  
  # Build data
  circuit = Circuit(*decl_val[1:])
  for stat in statements:
    #print list(stat)
    if stat[0] == import_:
      circuit.add_import(*stat[1:])
    elif stat[0] == gate:
      circuit.add_gate(*stat[1:])
    else:
      print stat
      raise Exception
  return circuit

def substitute(filename, args):
  # Parse for function declaration
  param_names = decl_stat.parseFile(filename)[2]
  params = {}
  assert len(param_names) == len(args), "Argument mismatch loading %s: len(%s) != len(%s)" % (filename, param_names, args)
  for name, val in zip(param_names, args):
    params[name] = val
  return process(filename, params)
