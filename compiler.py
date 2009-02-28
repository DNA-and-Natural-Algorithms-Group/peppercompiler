#!/usr/bin/env python
import pickle

from circuit_class import load_file, Circuit
from gate_class import Gate

def compiler(basename, args):
  """
  Start compiling a specification.
  
  Currently, it produces a .des file that can be used by sequence designers.
  """
  
  print "Compiling %s ..." % basename
  # Read in system (or component)
  system = load_file(basename, args)
  
  # TODO: Read in a forced sequences

  # Write the Zadeh-style design file
  outfile = file(basename+".des", "w")
  outfile.write("## Specification for %s compiled at: %s\n" % (basename, time.ctime()))
  system.output_nupack("", outfile)
  outfile.close()

  # Save compiler state to be reloaded when designer finishes
  save(system, basename+".save")
  print "System/component compiled into %s.des" % basename
  print "Run a designer on this and get an %s.mfe output" % basename
  print 'Finally run "python compiler.py %s.save" to finish compiling and run kinetics' % basename

def save(obj, filename):
  """Save an object for later finishing."""
  f = file(filename, "w")
  pickle.dump(obj, f)
  f.close()

def load(filename):
  """Load it back."""
  f = file(filename, "r")
  obj = pickle.load(f)
  f.close()
  return obj

if __name__ == "__main__":
  import sys
  import re
  import quickargs
  filename = sys.argv[1]
  args, keys = quickargs.get_args(sys.argv[2:])
  
  assert not keys, "Don't provide keywords to compiler.py"
  
  ## If we are starting compile specifying the entire filename
  p = re.match(r"(.*)\.(sys|comp)\Z", filename)
  if p:
    compiler(p.group(1), args)
  
  else: # If we are starting a compile with just the basename
    compiler(filename, args)
