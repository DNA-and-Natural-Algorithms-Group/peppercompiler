"""Nucleic Acid template design grammar"""

from template_class import Gate
from pyparsing import *
from var_substitute import process

## Some globals
# Pyparsing shortcuts
K = CaselessKeyword
S = Suppress
O = Optional
List = lambda x: Group(ZeroOrMore(x))
Map = lambda func: (lambda s,l,t: map(func, t) )
import string
lowers = string.lowercase
# syntax names and sets
# Codes for different subsets of nucleotides
NAcodes = "ACGTUNSWRYMKVHBD"
#func = "function"
func = "declare"
seq = "sequence"
sup_seq = "sequence"
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

seq_var = Group(Word(lowers, alphas+"_") + Optional("*", default=""))
seq_list = List(seq_var)

# Strand definition could be:
#  1) some basic sequence constraints like 5N or 2S or
#  2) a sequence variable (possibly complimented with *)
strand_const = seq_const | seq_var
strand_const_list = List(strand_const)

strand_var = var.copy()
strand_list = List(strand_var + S(Optional("+")))

struct_var = var.copy()
struct_list = List(struct_var + S(Optional("+")))

secondary_struct = Word( nums+"UH()+ " ) # I don't need to break it up

### TODO: allow ins and outs to be wc complements (i.e. seq_vars not just vars)
# function <gate> = <func name>(<params>): <inputs> -> <outputs>
params = O(S("(") + Group(delimitedList(var)) + S(")"), default=[])
func_stat = S(K(func)) + var + params + S(":") + var_list + S("->") + out_list

# sequence <name> = <constraints> : <length>
seq_stat  = S(K(seq))  + var + S("=") + seq_const_list + S(":") + integer

# sup-sequence <name> = <constraints / sequences> : <length>
sup_seq_stat = S(K(sup_seq)) + var + S("=") + strand_const_list + S(":") + integer

# strand <name> = <constraints / sequences> : <length>
strand_stat  = S(K(strand))  + var + S("=") + strand_const_list + S(":") + integer

# structure <optinoal mfe param> <name> = <strands> : <secondary structure>
mfe_info = O(  K("--no-mfe").setParseAction(lambda s,t,l: False) | \
               (S("--mfe=") + float_)  , default=1.0)
struct_stat = S(K(struct)) + mfe_info + var + S("=") + strand_list + S(":") + secondary_struct

# kin <inputs> -> <outputs>
kin_stat = S(K(kin)) + struct_list + S("->") + struct_list

statement = func_stat | seq_stat | sup_seq_stat | strand_stat | struct_stat | kin_stat

document = StringStart() + delimitedList(O(statement), "\n") + StringEnd()
document.ignore(pythonStyleComment)

def load_template(basename, args):
  ## replace args
  doc = substitute(basename, args)
  global gate # DEBUG
  gate = Gate()

  # Set parse actions to build data
  func_stat.setParseAction(lambda s,t,l: gate.add_function(l))
  seq_stat.setParseAction(lambda s,t,l: gate.add_sequence(l))
  sup_seq_stat.setParseAction(lambda s,t,l: gate.add_super_sequence(l))
  strand_stat.setParseAction(lambda s,t,l: gate.add_strand(l))
  struct_stat.setParseAction(lambda s,t,l: gate.add_structure(l))
  kin_stat.setParseAction(lambda s,t,l: gate.add_kinetics(l))
  
  seq_var.setParseAction(lambda s,t,l: (gate.seqs[l[0][0]] if l[0][1] == "" else ~gate.seqs[l[0][0]]))
  strand_var.setParseAction(lambda s,t,l: gate.strands[l[0]])
  struct_var.setParseAction(lambda s,t,l: gate.structs[l[0]])
  
  # Build data
  document.parseString(doc)
  assert gate.def_func
  return gate

def substitute(basename, args):
  filename = basename+".template"
  # Clear parse actions
  func_stat.setParseAction(lambda s,t,l: None)
  # Parse for function declaration
  print filename
  param_names = func_stat.parseFile(filename)[1]
  params = {}
  assert len(param_names) == len(args)
  for name, val in zip(param_names, args):
    params[name] = val
  return process(params, filename)






