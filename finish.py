#!/usr/bin/env python
from __future__ import division

import string
import math

from compiler import load
from kinetics import read_design, test_kinetics

from circuit_class import load_file, Circuit
from gate_class import Gate
from DNA_classes import wc

def finish(basename, **keys):
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
  seqs = read_design(basename+".mfe")
  
  # Apply designer results to our system
  apply_design(system, seqs)
  
  # Document all sequences, super-sequences, strands, and structures
  print "Writing sequences file: %s.seqs" % basename
  f = open(basename + ".seqs", "w")
  
  f.write("# Sequences\n")
  for name, seq in system.nupack_seqs.items():
    f.write("sequence %s\t%s\n" % (name, seq.seq))
  f.write("# Super-Sequences\n")
  for name, seq in system.sup_seqs.items():
    f.write("super-sequence %s\t%s\n" % (name, seq.seq))
  f.write("# Strands\n")
  for name, strand in system.strands.items():
    f.write("strand %s\t%s\n" % (name, strand.seq))
  f.write("# Structures\n")
  for name, struct in system.structs.items():
    f.write("struct %s\t%s\n" % (name, struct.seq))
  f.close()
  
  # Document all strands that will be used in the final experiment
  print 'Writing "stands to order" file: %s.strands' % basename
  f = open(basename + ".strands", "w")
  for name, strand in system.strands.items():
    if not strand.dummy:
      f.write("strand %s\t%s\n" % (name, strand.seq) )
  f.close()
  
  # TODO: Write a thermodynamic scorecard.
  
  # Run Kinetic tests
  print "Testing Kinetics with parameters: ", keys
  kinetic_rec(system, "", **keys)


def apply_design(system, seqs):
  """Assigns designed sequences and provided mfe structures to the respective objects."""
  # Assign all the designed sequences
  for name, seq in system.nupack_seqs.items():
    seq.seq  = seqs[name]
    assert len(seq.seq) == seq.length
    seq.wc.seq = wc(seq.seq)
    assert seq.wc.seq == seqs[name + "*"]
  
  for sup_seq in system.sup_seqs.values():
    sup_seq.seq  = string.join([seq.seq for seq in sup_seq.nupack_seqs], "")
    sup_seq.wc.seq = string.join([seq.seq for seq in sup_seq.wc.nupack_seqs], "")
  
  for strand in system.strands.values():
    strand.seq = string.join([seq.seq for seq in strand.nupack_seqs], "")
  
  for name, struct in system.structs.items():
    struct.seq = string.join([strand.seq for strand in struct.strands], "+")
    assert struct.seq == seqs[name], "Design is inconsistant! %s != %s" % (struct.seq, seqs[name])


# TODO: get rid of need for gate, so get rid of recurssion.
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
    try:
      num_trials, (coll_rate, coll_var), (for_rate, for_var, for_num, for_mean_time), (rev_rate, rev_var, rev_num, rev_mean_time) \
        = test_kinetics(kin, gate, **keys)
      # Process results
      coll_stddev = math.sqrt(coll_var) if coll_var != None else -1
      for_stddev = math.sqrt(for_var) if for_var != None else -1
      print "  Simulated %d trajectories, %.2f%% went forward." % (num_trials, 100*for_num/num_trials)
      print "  Collision Reaction Rate: %f (std-dev %f) (/M/s)" % (coll_rate, coll_stddev)
      print "  Forward Trajectory Rate: %f (std-dev %f) (/s) [Mean time: %f]" % (for_rate, for_stddev, for_mean_time)
    except:
      print "Something funny happened ... email Shawn <sligocki@gmail.com>"

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

if __name__ == "__main__":
  import sys
  import re
  import quickargs
  filename = sys.argv[1]
  args, keys = quickargs.get_args(sys.argv[2:])
  
  p = re.match(r"(.*)\.(save|mfe)\Z", filename)
  if p:
    basename = p.group(1)  # Circuit.save -> Circuit   or   Circuit.mfe -> Circuit
  else:
    basename = filename
  
  finish(basename, *args, **keys)

