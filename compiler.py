#!/usr/bin/env python
from __future__ import division

import sys
import pickle
import re
import time
import string

from circuit_class import load_file, Circuit
from gate_class import Gate
from kinetics import read_nupack, test_kinetics
import myStat as stat

def compiler(basename, args):
  """
  Start compiling a specification.
  
  Currently, it produces a .des file that can be used by sequence designers.
  """
  
  print "Compiling %s ..." % basename
  # Read in system (or component)
  system = load_file(basename, args)

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

def finish(basename, trials=24, num_proc=4, time=100000):
  """
  Finish compiling a specification.
  
  Currently, it produces .strands file of "strands to order"
                and tests the specified kinetic paths
  
  In the future it might test for bad kinetics, etc.
  """
  
  print "Finishing compilation of %s ..." % basename
  # Re-load the design system/component
  system = load(basename+".save")
  # Read results of DNA designer
  seqs, mfe_structs = read_nupack(basename+".mfe")
  
  # Apply designer results to our system
  apply_design_rec(system, "", seqs, mfe_structs)
  
  # Get all strands that will be used in the final experiment
  print 'Writing "stands to order" file: %s.strands' % basename
  strands = get_strands_rec(system)
  f = open(basename + ".strands", "w")
  for strand_name, strand_seq in strands:
    f.write("#Strand: %s\n%s\n\n" % (strand_name, strand_seq) )
  f.close()
  
  # Run Kinetic tests
  print "Testing Kinetics"
  print "Running %d trials across %d processes with max_time of %d." % (trials, num_proc, time)
  kinetic_rec(system, "", trials=trials, num_proc=num_proc, time=time)


def apply_design(gate, prefix, seqs, mfe_structs):
  """Assigns designed sequences and provided mfe structures to the respective objects."""
  # Assign all the designed sequences
  for seq in gate.reg_seqs.values() + gate.junk_seqs.values():
    seq.seq  = seqs[prefix + seq.name]
    assert len(seq.seq) == seq.length
    seq.wc.seq = seqs[prefix + seq.wc.name]
    assert len(seq.wc.seq) == seq.wc.length
  for sup_seq in gate.sup_seqs.values():
    sup_seq.seq  = string.join([seq.seq for seq in sup_seq.nupack_seqs], "")
    ~sup_seq #TODO: remove
    sup_seq.wc.seq = string.join([seq.seq for seq in sup_seq.wc.nupack_seqs], "")
  for strand in gate.strands.values():
    strand.seq = string.join([seq.seq for seq in strand.nupack_seqs], "")
  for struct in gate.structs.values():
    struct.seq = string.join([strand.seq for strand in struct.strands], "+")
    assert struct.seq == seqs[prefix + struct.name], "Design is inconsistant! %s != %s" % (struct.seq, seqs[prefix + struct.name])
  
  # Assign all the resulting mfe structures
  for struct in gate.structs.values():
    # TODO: Or should we just get it ourselves?  struct, dG = DNAfold(seq, temp)
    struct.mfe_struct = mfe_structs[prefix + struct.name]

def apply_design_rec(obj, prefix, seqs, mfe_structs):
  """Applies the results of a design to a system (which might be a single gate)."""
  # If it's actually a circuit, recurse.
  if isinstance(obj, Circuit):
    for gate_name, gate in obj.gates.items():
      apply_design_rec(gate, prefix + gate_name + "-", seqs, mfe_structs)
  
  # If it's a gate, use the gate code.
  elif isinstance(obj, Gate):
    return apply_design(obj, prefix, seqs, mfe_structs)
  else:
    raise Exception, 'Object "%r" is niether a Circuit or Gate.' % obj


def get_strands(gate):
  """Get all the strands that will be used in the final product."""
  # Collect the strand sequences 
  strands = [(strand.name, strand.seq) for strand in gate.strands.values() if not strand.dummy]
  return strands

def get_strands_rec(obj):
  """Get all strands in a system (which might be a single gate)."""
  # If it's actually a circuit, recurse.
  if isinstance(obj, Circuit):
    strands = []
    for gate_name, gate in obj.gates.items():
      strands += get_strands_rec(gate)
    return strands
  
  # If it's a gate, use the gate code.
  elif isinstance(obj, Gate):
    return get_strands(obj)
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
    frac, times, res = test_kinetics(kin, gate, **keys)
    # Process results
    print "%d%% finished. Mean time = %.1f  Std Dev = %.1f  Skewness = %.4f  Kurtosis = %.4f" \
          % (100*frac, stat.mean(times), stat.stddev(times), stat.skewness(times), stat.kurtosis(times))
    print "MLE Gamma distribution: k = %f, theta = %f" % stat.gamma_mle(times)

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
