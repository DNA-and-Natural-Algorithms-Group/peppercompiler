#!/usr/bin/env python
from __future__ import division

import sys
import pickle
import re
import time

from circuit_class import load_file, Circuit
from gate_class import Gate
from kinetics import read_nupack, test_kinetics
import myStat as stat

def compiler(basename, args):
  print "Compiling %s ..." % basename
  # Read in design (system or component)
  design = load_file(basename, args)

  # Write the Zadeh-style design file
  outfile = file(basename+".des", "w")
  outfile.write("## Specification for %s compiled at: %s\n" % (basename, time.ctime()))
  design.output_nupack("", outfile)
  outfile.close()

  # Save compiler state to be reloaded when designer finishes
  save(design, basename+".save")
  print "System/component compiled into %s.des" % basename
  print "Run a designer on this and get an %s.mfe output" % basename
  print 'Finally run "python compiler.py %s.save" to finish compiling and run kinetics' % basename

def finish(basename, trials=24, num_proc=4, time=100000):
  print "Finishing compilation of %s ..." % basename
  # Re-load the design system/component
  design = load(basename+".save")
  # Read results of DNA designer
  seqs, mfe_structs = read_nupack(basename+".mfe")
  
  # Get all strands that will be used in the final experiment
  print 'Writing "stands to order" file: %s.strands' % basename
  strands = get_strands_rec(design, "", seqs)
  f = open(basename + ".strands", "w")
  for strand_name, strand_seq in strands:
    f.write("#Strand: %s\n%s\n\n" % (strand_name, strand_seq) )
  f.close()
  
  # Run Kinetic tests
  print "Testing Kinetics"
  print "Running %d trials across %d processes with max_time of %d." % (trials, num_proc, time)
  kinetic_rec(design, "", seqs=seqs, mfe_structs=mfe_structs, trials=trials, num_proc=num_proc, time=time)


def get_strands(gate, prefix, seqs):
  """Get all the strands that will be used in the final product."""
  strands = []
  for strand in gate.strands.values():
    if not strand.dummy: # Ignore dummy strands
      strand_seq = ""
      for seq in strand.nupack_seqs:
        strand_seq += seqs[prefix + seq.name]
      strand_name = prefix + strand.name
      
      strands.append( (strand_name, strand_seq) )
  return strands

def get_strands_rec(obj, prefix, seqs):
  """Get all strands in a system (which might be a single gate)."""
  # If it's actually a circuit, recurse.
  if isinstance(obj, Circuit):
    strands = []
    for gate_name, gate in obj.gates.items():
      strands += get_strands_rec(gate, prefix + gate_name + "-", seqs)
    return strands
  
  # If it's a gate, use the gate code.
  elif isinstance(obj, Gate):
    return get_strands(obj, prefix, seqs)
  else:
    raise Exception, 'Object "%r" is niether a Circuit or Gate.' % obj


def kinetic(gate, prefix, **keys):
  """Test all kinetic pathways in a gate."""
  for kin in gate.kinetics.values():
    # Print testing info
    print
    print "kinetic", prefix[:-1], ":",
    for compl in kin.inputs:
      print compl.name,
    print "->",
    for compl in kin.outputs:
      print compl.name,
    print
    sys.stdout.flush()
    
    # Call Multistrand instances
    frac, times, res = test_kinetics(prefix, kin, **keys)
    # Process results
    # Currently: assumes an exponential distribuiton and finds the 
    #            99% confidence interval for the mean.
    alpha = 0.01
    low, mid, high = stat.exp_mean_interval(times, alpha)
    #print times
    print "%d%% finished. Mean time = %.0f (%.2f%% Confidence interval: %.0f < mean < %.0f)." \
          % (100*frac, mid, 100*(1-alpha), low, high)

def kinetic_rec(obj, prefix, **keys):
  """Run kinetic tests on gate (which might actually be a sub-circuit)."""
  # If it's actually a circuit, recurse.
  if isinstance(obj, Circuit):
    for gate_name, gate in obj.gates.items():
      kinetic_rec(gate, prefix + gate_name + "-", **keys)
  
  # If it's a gate, use the gate code.
  elif isinstance(obj, Gate):
    kinetic(obj, prefix, **keys)
  else:
    raise Exception, 'Object "%r" is niether a Circuit or Gate.' % obj


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
