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
H = Hidden = lambda x: Empty().setParseAction(lambda s,t,l: x)  # A hidden field, it tags an entry
List = lambda x: Group(ZeroOrMore(x))
Map = lambda func: (lambda s,l,t: map(func, t) )

import string
lowers = string.lowercase

# Syntax names and sets
# Codes for different subsets of nucleotides
NAcodes = "ACGTUNSWRYMKVHBD"

decl = "declare"
comp = "component"
seq = "sequence"
sup_seq = "sequence"; sup_seq_key = "sup-sequence"
strand = "strand"
struct = "structure"
kin = "kinetic"

# Don't ignore newlines!
ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
integer = Word(nums).setParseAction(Map(int))
float_ = Word(nums+"+-.eE").setParseAction(Map(float))

# Signals used in the declare line seq
sig_list = List(var + S(O("+")))

# Sequence const could be ?N, 3N or N
seq_const = Group(( "?" | Optional(integer, default=1) ) + Word(NAcodes, exact=1))
seq_const_list = List(seq_const)

seq_name = Word(lowers, alphanums+"_") # Sequence name starts with lower case
seq_var = Group(H("Sequence") + Group(seq_name + Optional("*", default="")))
seq_list = List(seq_var)

# Strand definition could be:
#  1) Some basic sequence constraints like 5N or 2S or
#  2) A sequence variable (possibly complimented with *) (Must start with lowercase letter)
strand_const = seq_const | seq_var
strand_const_list = List(strand_const)

strand_var = var;  strand_list = List(strand_var + S(Optional("+")))
struct_var = var;  struct_list = List(struct_var + S(Optional("+")))

# Secondary structure can be dot-paren, extended dot-paren or HU notation.
## TODO: deal with errors. Super-confusing error messages right now.
exDotParen = Word(nums + ".()+ " ).setParseAction(Map(extended2dotParen))
HUnotation = Word(nums + "UH()+ ").setParseAction(Map(HU2dotParen))
secondary_struct = exDotParen | HUnotation


### TODO: allow ins and outs to be wc complements (i.e. seq_vars not just vars)
# declare component <gate name>(<params>): <inputs> -> <outputs>
params = O(S("(") + Group(delimitedList(var)) + S(")"), default=[])
decl_stat = K(decl) + S(comp) + var + params + S(":") + sig_list + S("->") + sig_list
# sequence <name> = <constraints> : <length>
seq_stat  = K(seq)  + seq_name + S("=") + seq_const_list + S(":") + integer
# sup-sequence <name> = <constraints / sequences> : <length>
sup_seq_stat = K(sup_seq).setParseAction(lambda s,t,l: sup_seq_key) + \
               seq_name + S("=") + strand_const_list + S(":") + integer
# strand <name> = <constraints / sequences> : <length>
strand_stat  = K(strand) + O("[dummy]", default="") + strand_var + S("=") + strand_const_list + S(":") + integer
# structure <optinoal opt param> <name> = <strands> : <secondary structure>
opt = Optional(   K("[no-opt]").setParseAction(lambda s,t,l: False) | \
                ( Suppress("[") + float_ + Suppress("nt]") ),
                  default=1.0)
struct_stat = K(struct) + opt + struct_var + S("=") + strand_list + S(":") + secondary_struct
# kin <inputs> -> <outputs>
kin_stat = K(kin) + struct_list + S("->") + struct_list

statement = seq_stat | sup_seq_stat | strand_stat | struct_stat | kin_stat
document = StringStart() + ZeroOrMore(S("\n")) + Group(decl_stat) + S("\n") + \
           Group(delimitedList(O(Group(statement)), "\n")) + StringEnd()
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
    decl_val, statements = document.parseString(doc)
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
