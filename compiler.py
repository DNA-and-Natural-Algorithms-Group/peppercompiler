#!/usr/bin/env python
from __future__ import division

import sys
import pickle
import re
import time
import myStat as stat

from circuit_class import load_gate, Circuit
from template_class import Gate
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
  obj = load(basename+".save")
  # Read results
  seqs, mfe_structs = read_nupack(basename+".mfe")
  
  # Helper functions
  def kin(obj, prefix=""):
    """Run kinetic tests on object (could be ciruit or gate)."""
    if isinstance(obj, Circuit):
      return kin_circuit(obj, prefix)
    elif isinstance(obj, Gate):
      return kin_gate(obj, prefix)
    else:
      raise Exception

  def kin_circuit(circuit, prefix):
    """Run kinetics on all gates (and subcircuits) a given circuit."""
    for obj_name, obj in circuit.gates.items():
      kin(obj, prefix + obj_name + "-")

  def kin_gate(gate, prefix):
    """Test all kinetic pathways in a gate."""
    for kin in gate.kinetics.values():
      # Print testing info
      print
      print "kinetic", prefix, ":",
      for compl in kin.inputs:
        print compl.name,
      print "->",
      for compl in kin.outputs:
        print compl.name,
      print
      sys.stdout.flush()
      # Call Multistrand instances
      frac, times, res = test_kinetics(prefix, kin, seqs, mfe_structs, trials=trials, num_proc=num_proc, time=time)
      # TODO: process results
      alpha = 0.01
      low, mid, high = stat.exp_mean_interval(times, alpha)
      print times
      print "%d%% finished. Mean time = %.0f (%.2f%% Confidence interval: %.0f < mean < %.0f)." \
            % (100*frac, mid, 100*(1-alpha), low, high)
  # End of helper functions
  
  # Prepare for Schaeffer's Multistrand
  kin(obj)

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
