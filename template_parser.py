"""Nucleic Acid template design grammar"""

from template_class import Gate
from pyparsing import *
from var_substitute import process

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
func = "declare"
seq = "sequence"
sup_seq = "sequence"; sup_seq_key = "sup-sequence"
strand = "strand"
struct = "structure"
kin = "kin"

ParserElement.setDefaultWhitespaceChars(" \t")


## Define Grammar
var = Word(alphas, alphanums+"_") # Variable name
var_list = List(var + S(O("+")))  # Space sep list of variable names optional plus
integer = Word(nums).setParseAction(Map(int))
float_ = Word(nums+"+-.eE").setParseAction(Map(float))

out = var + S("(" + var + ")")
out_list = List(out + S(O("+")))

# Sequence const could be ?N, 3N or N
seq_const = Group(( "?" | Optional(integer, default=1) ) + Word(NAcodes))
seq_const_list = List(seq_const)

seq_var = Group(H("Sequence") + Group(Word(lowers, alphas+"_") + Optional("*", default="")))
seq_list = List(seq_var)

# Strand definition could be:
#  1) Some basic sequence constraints like 5N or 2S or
#  2) A sequence variable (possibly complimented with *) (Must start with lowercase letter)
strand_const = seq_const | seq_var
strand_const_list = List(strand_const)

strand_var = var
strand_list = List(strand_var + S(Optional("+")))

struct_var = var
struct_list = List(struct_var + S(Optional("+")))

secondary_struct = Word( nums+"UH()+ " ) # I don't need to break it up

### TODO: allow ins and outs to be wc complements (i.e. seq_vars not just vars)
# function <gate> = <func name>(<params>): <inputs> -> <outputs>
params = O(S("(") + Group(delimitedList(var)) + S(")"), default=[])
func_stat = K(func) + var + params + S(":") + var_list + S("->") + out_list

# sequence <name> = <constraints> : <length>
seq_stat  = K(seq)  + var + S("=") + seq_const_list + S(":") + integer

# sup-sequence <name> = <constraints / sequences> : <length>
sup_seq_stat = K(sup_seq).setParseAction(lambda s,t,l: sup_seq_key) + var + S("=") + strand_const_list + S(":") + integer

# strand <name> = <constraints / sequences> : <length>
strand_stat  = K(strand)  + var + S("=") + strand_const_list + S(":") + integer

# structure <optinoal mfe param> <name> = <strands> : <secondary structure>
mfe_info = O(  K("--no-mfe").setParseAction(lambda s,t,l: False) | \
               (S("--mfe=") + float_)  , default=1.0)
struct_stat = K(struct) + mfe_info + var + S("=") + strand_list + S(":") + secondary_struct

# kin <inputs> -> <outputs>
kin_stat = K(kin) + struct_list + S("->") + struct_list

statement = func_stat | seq_stat | sup_seq_stat | strand_stat | struct_stat | kin_stat

document = StringStart() + delimitedList(O(Group(statement)), "\n") + StringEnd()
document.ignore(pythonStyleComment)

def load_template(basename, args):
  ## replace args
  doc = substitute(basename, args)
  global gate # DEBUG
  gate = Gate()

  # Build data
  statements = document.parseString(doc)
  for stat in statements:
    #print list(stat)
    if stat[0] == func:
      gate.add_function(stat[1:])
    elif stat[0] == seq:
      gate.add_sequence(stat[1:])
    elif stat[0] == sup_seq_key:
      gate.add_super_sequence(stat[1:])
    elif stat[0] == strand:
      gate.add_strand(stat[1:])
    elif stat[0] == struct:
      gate.add_structure(stat[1:])
    elif stat[0] == kin:
      gate.add_kinetics(stat[1:])
    else:
      print stat
      raise Exception
  assert gate.def_func
  return gate

def substitute(basename, args):
  filename = basename+".template"
  # Parse for function declaration
  #print filename
  param_names = func_stat.parseFile(filename)[2]
  params = {}
  #print param_names, args
  assert len(param_names) == len(args)
  for name, val in zip(param_names, args):
    params[name] = val
  return process(params, filename)






