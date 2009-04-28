#!/usr/bin/env python
import re
import os
import pickle
import time

from circuit_class import load_file, Circuit
from gate_class import Gate
from utils import match

def parse_fixed(line):
  """Parse a line in the fixed file."""
  type_, name, seq = match(r"(\w+) ([\w_-]+) ([ATCG]+)(?: #.*)?", line)
  return type_, name, seq

def load_fixed(filename):
  """Load a file of sequences to fix."""
  if os.path.isfile(filename):
    f = open(filename, "r")
    return [parse_fixed(line) for line in f if not re.match(r"\s*(#.*)?\s*\Z", line)]
  else:
    return []

def compiler(basename, args, outputname, savename, fixed_file=None):
  """
  Start compiling a specification.
  
  Currently, it produces a .des file that can be used by sequence designers.
  """
  
  print "Compiling '%s' ..." % basename
  # Read in system (or component)
  system = load_file(basename, args)
  
  if fixed_file:
    print "Fixing sequences from file '%s'" % fixed_file
    fixed_sequences = load_fixed(fixed_file)
    for type_, name, fixed_seq in fixed_sequences:
      if type_ in ("sequence", "super-sequence"):
        system.seqs[name].fix_seq( fixed_seq )
      elif type_ == "strand":
        system.strands[name].fix_seq( fixed_seq )
      elif type_ == "structure":
        system.structs[name].fix_seq( fixed_seq )
  

  # Write the Zadeh-style design file
  print "System/component compiled into '%s'" % outputname
  outfile = open(outputname, "w")
  outfile.write("## Specification for %s compiled at: %s\n" % (basename, time.ctime()))
  system.output_nupack("", outfile)
  outfile.close()

  # Save compiler state to be reloaded when designer finishes
  print "Compiler state saved into '%s'" % savename
  save(system, savename)
  print "Run a designer on '%s' and process the result with python finish.py" % outputname

def save(obj, filename):
  """Save an object for later finishing."""
  f = open(filename, "w")
  pickle.dump(obj, f)
  f.close()

def load(filename):
  """Load it back."""
  f = open(filename, "r")
  obj = pickle.load(f)
  f.close()
  return obj

if __name__ == "__main__":
  from optparse import OptionParser
  
  # Parse command line options.
  usage = "usage: %prog [options] BASENAME [parameters ...]"
  parser = OptionParser(usage=usage)
  #parser.set_defaults(verbose=True)
  #parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
  # TODO: implement quiet
  #parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
  parser.add_option("--fixed", help="Fix specific sequences listed in FILE", metavar="FILE")
  parser.add_option("--output", help="Output file [defaults to BASENAME.des]", metavar="FILE")
  parser.add_option("--save", help="Saved state file [defaults to BASENAME.save]", metavar="FILE")
  (options, args) = parser.parse_args()
  
  # Get basename of input specification
  if len(args) < 1:
    parser.error("missing required argument BASENAME")
  basename = args[0]
  # Infer the basename if a full filename is given
  p = re.match(r"(.*)\.(sys|comp)\Z", basename)
  if p:
    basename = p.group(1)
  
  # Set filename defaults
  if not options.output:
    options.output = basename + ".des"
  if not options.save:
    options.save = basename + ".save"
  
  # Eval remaining arguments
  import quickargs
  args, keys = quickargs.get_args(args[1:])
  assert not keys, "Don't provide keywords to compiler.py"
  
  compiler(basename, args, options.output, options.save, options.fixed)
