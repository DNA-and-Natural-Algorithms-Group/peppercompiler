"""Nucleic Acid template design grammar"""

from nupack_in_class import Spec
from pyparsing import *

## Some globals
# Pyparsing shortcuts
K = CaselessKeyword
S = Suppress
O = Optional
List = lambda x: Group(OneOrMore(x))
Map = lambda func: (lambda s,l,t: map(func, t) )
# syntax names and sets
# Codes for different subsets of nucleotides
NAcodes = "ACGTUNSWRYMKVHBD"
struct = "structure"
seq = "sequence"

ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas+"_-", alphanums+"_-") # Variable name

integer = Word(nums).setParseAction(Map(int))
float_ = Word(nums+"-.").setParseAction(Map(float))

# Sequence const could be ?N, 3N or N
seq_const = Group(Optional(integer, default=1) + Word(NAcodes))
seq_const_list = List(seq_const)

seq_var = Group(var+Optional("*", default=""))
seq_list = List(seq_var)

struct_var = var.copy()
struct_list = List(struct_var)

### TODO: Break it up
secondary_struct = Word( nums+"UH()+ " )

# structure <name> = <secondary structure>
struct_stat = S(K(struct)) + var + S("=") + secondary_struct
# sequence <name> = <constraints>
seq_stat  = S(K(seq))  + var + S("=") + seq_const_list
# prevent <list of structs> < <float>
prevent_stat = S(K("prevent")) + struct_list + S("<") + float_
# <struct name> : <list of seqs>
apply_stat = struct_var + S(":") + seq_list
# <struct name> < <float>
obj_stat = struct_var + S("<") + float_

statement = struct_stat | seq_stat | prevent_stat | apply_stat | obj_stat

document = StringStart() + delimitedList(O(statement), "\n") + StringEnd()
document.ignore(pythonStyleComment)

def ifelse(cond, true_part, false_part):
  if cond:
    return true_part()
  else:
    return false_part()
                                  

def load_design(filename):
  spec = Spec()

  # Set parse actions to build data
  struct_stat.setParseAction(lambda s,t,l: spec.add_structure(l))
  seq_stat.setParseAction(lambda s,t,l: spec.add_sequence(l))
  apply_stat.setParseAction(lambda s,t,l: spec.add_apply(l))
  # ignore the prevent and objective statements ...
  seq_var.setParseAction(lambda s,t,l: ifelse(l[0][1] == "", lambda: spec.seqs[l[0][0]], lambda: ~spec.seqs[l[0][0]]))
  # Python 2.5 is ok with: seq_var.setParseAction(lambda s,t,l: (spec.seqs[l[0][0]] if l[0][1] == "" else ~spec.seqs[l[0][0]]))
  struct_var.setParseAction(lambda s,t,l: spec.structs[l[0]])
  
  # Build data
  document.parseFile(filename)
  return spec


