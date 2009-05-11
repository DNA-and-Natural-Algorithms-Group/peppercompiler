import sys

from circuit_class import Circuit
from var_substitute import process

from pyparsing import *

## Some globals
K = CaselessKeyword
S = Suppress
O = Optional
Map = lambda func: (lambda s, l, t: map(func, t) )  # A useful mapping function

def List(expr, delim=""):
  """My delimited list. Allows for length zero list and uses no delimiter by default."""
  if not delim:
    return Group(ZeroOrMore(expr))
  else:
    return Group(Optional(expr + ZeroOrMore( Suppress(delim) + expr)))

def Flag(expr):
  """A flag identifier. It is either present or not and returns True or False."""
  p = Optional(expr)
  p.setParseAction(lambda s, l, t: bool(t))
  return p

decl = "declare"
system = "system"
import_ = "import"
gate = "component"

# Don't ignore newlines!
ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
signal = Group(var + Flag("*"))
signal_list = List(signal, "+")

path = Word(alphanums+".-_/~") # Path name in a directory structure
py_chars = printables.replace(",", "").replace(")", "")
python_object = Word(py_chars, py_chars+" ").setParseAction(Map(eval))


# declare system <cicuit name>(<params>): <inputs> -> <outputs>
decl_params = O(S("(") + List(var, ",") + S(")"), default=[])
decl_stat = K(decl) + S(system) + var + decl_params + S(":") + signal_list + S("->") + signal_list

# import Adder, HalfAdder5 as HalfAdder, templates/Crossing_Gates/LastAdder
import_stat = K(import_) + delimitedList(Group(path + O(S("as") + var, default=None)))

# gate <name> = <template name>(<params>): <inputs> -> <outputs>
gate_params = O( S("(") + List(python_object, ",") + S(")") , default=[])
gate_stat = K(gate) + var + S("=") + var + gate_params + S(":") + signal_list + S("->") + signal_list


statement = import_stat | gate_stat

document = StringStart() + Group(decl_stat) + S("\n") + \
           List(O(Group(statement)), delim="\n") + StringEnd()
document.ignore(pythonStyleComment)



def load_circuit(filename, args, path):
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
    declare, statements = document.parseString(doc, parseAll=True)
  except ParseBaseException, e:
    print
    print doc
    print "Parsing error in circuit:", filename
    print e
    sys.exit(1)
  
  x, name, params, inputs, outputs = declare
  # Build data
  circuit = Circuit(path, name, params)
  for stat in statements:
    #print list(stat)
    if stat[0] == import_:
      circuit.add_import(*stat[1:])
    elif stat[0] == gate:
      circuit.add_gate(*stat[1:])
    else:
      raise Exception, "Unexpected statement:\n" + stat
  circuit.add_IO(inputs, outputs)
  return circuit

def substitute(filename, args):
  # Parse for function declaration
  param_names = decl_stat.parseFile(filename)[2]
  params = {}
  assert len(param_names) == len(args), "Argument mismatch loading %s: len(%s) != len(%s)" % (filename, param_names, args)
  for name, val in zip(param_names, args):
    params[name] = val
  return process(filename, params)
