"""DNA component design grammar with regular expressions"""

import re
import sys

from HU2dotParen import extended2dotParen, HU2dotParen
from DNA_classes import group
from utils import error, match

#TODO: use these
sequence_flag = "sequence"
domains_flag = "domains"
nucleotide_flag = "nucleotide"

def parse_constraints(constraints):
  """Parse a list of constraints"""
  if not match(r'\A(("[\w?\s]+"|[\w-]+[*]?)(\Z|\s+))*\Z', constraints):
    error("Invalid sequence constraints format:\n"
          'Valid example: "5N" seq1 "SNNCTB" Toe_SEQ\n'
          "Was:           %s" % constraints)
  constraints = re.findall(r'("[\w?\s]+"|[\w-]+[*]?)', constraints)
  constraints = map(parse_constraint, constraints)
  return constraints

def parse_constraint(word):
  """Helper function for parsing sequence constraints"""
  m = match(r'"([\w\s?-]+)"', word)
  if m:
    constraint = m.group(1)
    constraint = re.sub(r"\s*", "", constraint)
    parts = re.split(r'(\d+|\?)', constraint)
    
    nts = []
    factor = 1 # Keep track of what the last factor was
    for part in parts:
      if not part:
        continue # Skip empty parts
      if re.match(r'\d+', part): # An integer is a factor for the next symbol
        factor = int(part)
      elif part == "?":
        factor = "?"
      else:
        nts.append( [factor, part[0]] ) # Factor of the first guy
        factor = 1
        nts += [ [1, x] for x in part[1:] ] # And one of the rest until we get to the next break
    return [nucleotide_flag, nts]
  
  m = match(r"([\w-]+)([*])?", word)
  if m:
    sequence_name, reverse = m.group(1, 2)
    reverse = bool(reverse)
    return [sequence_flag, [sequence_name, reverse]]

def parse_signal(signal):
  """Helper function for parsing input and output signals"""
  m = match(r"([\w-]+)(\*)?(\(([\w-]+)\))?", signal)
  if not m:
    error("Invalid signal format:\n"
          "Should be: <sequence name>   or   <sequence name>(<structure name>)\n"
          "Was:       %s" % signal)
  sequence, reverse, structure = m.group(1, 2, 4)
  reverse = bool(reverse)
  return [[sequence, reverse], structure]

def parse_declare_statement(line):
  """Parse declare statement"""
  m = match(r"declare component ([\w-]+)(\((.*)\))?:(.*) -> (.*)", line)
  if not m:
    error("Invalid component declare statement format:\n"
          "Should be: declare component <name>: <inputs> -> <outputs>\n"
          "or:        declare component <name>(<params>): <inputs> -> <outputs>\n"
          "Was:       %s" % line)
  name, params, inputs, outputs = m.group(1, 3, 4, 5)
  if params:
    params = params.split(",")
    params = [x.strip() for x in params if x.strip()]
  else:
    params = []
  inputs = inputs.split("+")
  inputs = [x.strip() for x in inputs if x.strip()]
  inputs = map(parse_signal, inputs)
  outputs = outputs.split("+")
  outputs = [x.strip() for x in outputs if x.strip()]
  outputs = map(parse_signal, outputs)
  return [name, params, inputs, outputs]

def parse_sequence_statement(line):
  """Parse sequence statements"""
  m = match(r"sequence ([\w-]+) = ([^:\s]+)( : (\d+))?", line)
  if not m:
    error("Invalid sequence statement format:\n"
          "Should be: sequence <name> = <constraints>\n"
          "or:        sequence <name> = <constraints> : <length>\n"
          "Was:       %s" % line)
  name, template, length = m.group(1, 2, 4)
  if not set(template).issubset( set(group.keys()) ):
    error("Sequence's constraint template must be in allowed alphabet (%r).\nLine: %s" % (group.keys(), line))
  if length:
    length = int(length)
  return [name, template, length]

# TODO: make parse_general_sequence_statement()
def parse_general_sequence_statement(line):
  """Parse sequence or super-sequence statement"""
  m = match(r"sequence ([\w-]+) = ([^:]+)( : (\d+))?", line)
  if not m:
    error("Invalid sequence statement format:\n"
          "Should be: sequence <name> = <constraints>\n"
          "or:        sequence <name> = <constraints> : <length>\n"
          "Was:       %s" % line)
  name, constraints, length = m.group(1, 2, 4)
  constraints = parse_constraints(constraints)
  if length:
    length = int(length)
  return [name, constraints, length]

# TODO: doesn't allow nucleotides
def parse_super_sequence_statement(line):
  """Parse super-sequence statements"""
  m = match(r"sup(er)?-sequence ([\w-]+) = ([^:]+)( : (\d+))?", line)
  if not m:
    error("Invalid super-sequence statement format:\n"
          "Should be: super-sequence <name> = <constraints>\n"
          "or:        super-sequence <name> = <constraints> : <length>\n"
          "Was:       %s" % line)
  name, seq_names, length = m.group(2, 3, 5)
  seq_names = seq_names.split()
  if length:
    length = int(length)
  return [name, seq_names, length]

def parse_strand_statement(line):
  """Parse strand statements"""
  m = match(r"strand (\[dummy\] )?([\w-]+) = ([^:]+)( : (\d+))?", line)
  if not m:
    error("Invalid strand statement format:\n"
          "Should be: strand <name> = <constraints>\n"
          "or:        strand <name> = <constraints> : <length>\n"
          "or:        strand [dummy] <name> = <constraints> : <length>\n"
          "Was:       %s" % line)
  dummy, name, constraints, length = m.group(1, 2, 3, 5)
  dummy = bool(dummy)
  constraints = parse_constraints(constraints)
  if length:
    length = int(length)
  return [dummy, name, constraints, length]

def parse_structure_statement(line):
  """Parse structure statements"""
  m = match(r"structure( \[([\wd\.]+)nt\])? ([\w-]+) = ([^:]+) :( domain)? (.+)", line)
  if not m:
    error("Invalid structure statement format:\n"
          "Should be: structure <name> = <strand names> : <secondary structure>\n"
          "or:        structure [<parameters>] <name> = <strand names> : <secondary structure>\n"
          "or:        structure [<parameters>] <name> = <strand names> : domain <secondary structure>\n"
          "Was:       %s" % line)
  opt, name, strand_names, domain, struct = m.group(2, 3, 4, 5, 6)
  if opt:
    opt = float(opt)
  else:
    opt = 1.0
  strand_names = strand_names.split("+")
  strand_names = [sname.strip() for sname in strand_names]  # Clean off whitespace
  domain = bool(domain)
  if "U" in struct or "H" in struct:
    struct = HU2dotParen(struct)
  else:
    struct = extended2dotParen(struct)
  return [opt, name, strand_names, [domain, struct]]

def parse_kinetic_statement(line):
  """Parse kinetic statements"""
  m = match(r"kinetic (.*) -> (.*)", line)
  if not m:
    error("Invalid strand statement format:\n"
          "Should be: kinetic <inputs> -> <outputs>\n"
          "Was:       %s" % line)
  inputs, outputs = m.group(1, 2)
  inputs = inputs.split("+")
  inputs = [x.strip() for x in inputs if x.strip()]
  outputs = outputs.split("+")
  outputs = [x.strip() for x in outputs if x.strip()]
  return [inputs, outputs]

def parse_equal_statement(line):
  """Parse super-sequence statements"""
  seq_names = line.split()[1:]
  return seq_names

