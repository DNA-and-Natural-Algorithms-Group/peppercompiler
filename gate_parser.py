"""Nucleic Acid gate design grammar"""

import sys

from HU2dotParen import extended2dotParen, HU2dotParen
from gate_class import Gate
from var_substitute import process

from pyparsing import *

## Some globals
# Pyparsing shortcuts
K = CaselessKeyword
S = Suppress
O = Optional
H = Hidden = lambda x: Empty().setParseAction(lambda s,l,t: x)  # A hidden field, it tags an entry
Map = lambda func: (lambda s,l,t: map(func, t) )

def List(expr, delim=""):
  """My delimited list. Allows for length zero list and uses no delimiter by default."""
  if not delim:
    return Group(ZeroOrMore(expr))
  else:
    return Group(Optional(expr + ZeroOrMore( Suppress(delim) + expr)))

def Flag(expr):
  """A flag identifier. It is either present or not and returns True or False."""
  p = Optional(expr)
  p.setParseAction(lambda s,l,t: bool(t))
  return p

import string
lowers = string.lowercase

# Syntax names and sets
# Codes for different subsets of nucleotides
NAcodes = "ACGTUNSWRYMKVHBD"

decl = "declare"
comp = "component"
seq = "sequence"
sup_seq = "sup-sequence"; sup_seq_key = "sup-sequence"
strand = "strand"
struct = "structure"
kin = "kinetic"

# Don't ignore newlines!
ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
integer = Word(nums).setParseAction(Map(int))
float_ = Word(nums+"+-.eE").setParseAction(Map(float))

# Sequence const could be ?N, 3N or N
seq_const = Group(( "?" | Optional(integer, default=1) ) + Word(NAcodes, exact=1))
seq_const_list = List(seq_const)

seq_name = Word(lowers, alphanums+"_") # Sequence name starts with lower case
seq_var = Group(H("Sequence") + Group(seq_name + Optional("*", default="")))

# Signals used in the declare line seq
signal_var = Group(seq_name + Optional("*", default=""))

# Strand definition could be:
#  1) Some basic sequence constraints like 5N or 2S or
#  2) A sequence variable (possibly complimented with *) (Must start with lowercase letter)
strand_const = seq_const | seq_var

strand_var = var
struct_var = var

# Secondary structure can be dot-paren, extended dot-paren or HU notation.
# 'domain' keyword means that the helix/unpaired segments are by domain
## TODO: deal with errors. Super-confusing error messages right now.
exDotParen = Word(nums + ".()+ " ).setParseAction(Map(extended2dotParen))
HUnotation = Word(nums + "UH()+ ").setParseAction(Map(HU2dotParen))
secondary_struct = Group( Flag("domain") + (exDotParen | HUnotation) )


# declare component <gate name>(<params>): <inputs> -> <outputs>
params = O(S("(") + List(var, ",") + S(")"), default=[])
decl_stat = K(decl) + S(comp) + var + params + S(":") + List(signal_var, O("+")) + S("->") + List(signal_var, O("+"))

# sequence <name> = <constraints> : <length>
seq_stat  = K(seq)  + seq_name + S("=") + List(seq_const) + O(S(":") + integer, default=None)

# sup-sequence <name> = <constraints / sequences> : <length>
sup_seq_stat = K(sup_seq) + \
               seq_name + S("=") + List(strand_const) + O(S(":") + integer, default=None)

# strand <name> = <constraints / sequences> : <length>
strand_stat  = K(strand) + Flag("[dummy]") + strand_var + S("=") + List(strand_const) + O(S(":") + integer, default=None)

# structure <optinoal opt param> <name> = <strands> : <secondary structure>
opt = Optional(   K("[no-opt]").setParseAction(lambda s,t,l: False) | \
                ( Suppress("[") + float_ + Suppress("nt]") ),
                  default=1.0)
struct_stat = K(struct) + opt + struct_var + S("=") + List(strand_var, O("+")) + S(":") + secondary_struct

# kin <inputs> -> <outputs>
kin_stat = K(kin) + List(struct_var, O("+")) + S("->") + List(struct_var, O("+"))


statement = seq_stat | sup_seq_stat | strand_stat | struct_stat | kin_stat
document = StringStart() + ZeroOrMore(S("\n")) + Group(decl_stat) + S("\n") + \
           List(O(Group(statement)), "\n") + StringEnd()
document.ignore(pythonStyleComment)



def load_gate(filename, args):
  """Load component file"""
  try:
    # Open file and do parameter substitution
    doc = substitute(filename, args)
  except ParseBaseException, e:
    print
    print "Parsing error in component:", filename
    print e
    sys.exit(1)
    
  try:
    # Load data
    decl_val, statements = document.parseString(doc, parseAll=True)
  except ParseBaseException, e:
    print
    print doc
    print "Parsing error in component:", filename
    print e
    sys.exit(1)
  
  # Build data
  gate = Gate(*decl_val[1:])
  for stat in statements:
    #print list(stat)
    if stat[0] == seq:
      gate.add_sequence(*stat[1:])
    elif stat[0] == sup_seq_key:
      gate.add_super_sequence(*stat[1:])
    elif stat[0] == strand:
      gate.add_strand(*stat[1:])
    elif stat[0] == struct:
      gate.add_structure(*stat[1:])
    elif stat[0] == kin:
      gate.add_kinetics(*stat[1:])
    else:
      print stat
      raise Exception
  return gate

def substitute(filename, args):
  # Parse for function declaration
  param_names = decl_stat.parseFile(filename)[2]
  params = {}
  assert len(param_names) == len(args), "Argument mismatch loading %s: len(%s) != len(%s)" % (filename, param_names, args)
  for name, val in zip(param_names, args):
    params[name] = val
  return process(filename, params)
