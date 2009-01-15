#!/usr/bin/env python
from __future__ import division

import sys
import pickle
import re
import time

from circuit_class import load_gate, Circuit
from kinetics import read_nupack, test_kinetics

def compiler(basename, args):
  print "Compiling %s ..." % basename
  # Read in circuit design
  circuit = load_gate(basename, args)
  # TODO: allow circuit to be a gate
  if not isinstance(circuit, Circuit):
    print "Warning: compiling Gates are not completely supported yet."

  # Write the Zadeh-style design file
  outfile = file(basename+".des", "w")
  outfile.write("## Specification for %s compiled at: %s\n" % (basename, time.ctime()))
  circuit.output_nupack("", outfile)
  outfile.close()

  # Save compiler state to be reloaded when designer finishes
  save(circuit, basename+".save")
  print "System/component compiled into %s.des" % basename
  print "Run a designer on this and get an %s.mfe output" % basename
  print 'Finally run "python compiler.py %s.save" to finish compiling and run kinetics' % basename

def finish(basename, trials=24, num_proc=4, time=100000):
  print "Finishing compilation of %s ..." % basename
  print "Running %d trials across %d processes with max_time of %d." % (trials, num_proc, time)
  circuit = load(basename+".save")
  # Read results
  seqs, mfe_structs = read_nupack(basename+".mfe")
  
  # Prepare for Schaeffer's Multistrand
  # TODO: allow circuit to be a gate
  if not isinstance(circuit, Circuit):
    print "Error: compiling Gates is not completely supported yet."
    return 1
  for gate_name, gate in circuit.gates.items():
    for kin in gate.kinetics.values():
      # Print testing info
      print "kinetic", gate_name, ":",
      for compl in kin.inputs:
        print compl.name,
      print "->",
      for compl in kin.outputs:
        print compl.name,
      print
      sys.stdout.flush()
      # Call Multistrand instances
      frac, times, res = test_kinetics(gate_name, kin, seqs, mfe_structs, trials=trials, num_proc=num+proc, time=time)
      if times:
        ave = sum(times) / len(times)
      else:
        ave = 0
      # TODO: process results
      print "Result:", frac, ave

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
  import quickargs
  filename = sys.argv[1]
  args, keys = quickargs.get_args(sys.argv[2:])
  
  ## If we are starting compile specifying the entire filename
  p = re.match(r"(.*)\.(sys|comp)\Z", filename)
  if p:
    compiler(p.group(1), args)
  
  else: # If we are finishing compile
    p = re.match(r"(.*)\.save\Z", filename)
    if p:
      finish(p.group(1), **keys) # Finish started process
    
    else: # If we are starting a compile with just the basename
      compiler(filename, args)
