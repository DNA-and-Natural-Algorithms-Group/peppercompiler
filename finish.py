#!/usr/bin/env python
from __future__ import division

import sys
import string

from compiler import load
from kinetics import read_design, test_kinetics

from circuit_class import Circuit
from gate_class import Gate
from DNA_classes import wc
import myStat as stat

def get_gates(obj, prefix=""):
  """Return all gates and the prefixes showing how to reach them."""
  if isinstance(obj, Gate):
    return [(obj, prefix)]
  elif isinstance(obj, Circuit):
    gates = []
    for gate_name, gate in obj.gates.items():
      gates += get_gates(gate, prefix + gate_name + "-")
    return gates
  else:
    raise Exception, 'Object "%r" is neither a Circuit or Gate.' % obj

def finish(savename, designname, seqsname, strandsname, run_kin, cleanup, keys):
  """
  Finish compiling a specification.
  
  Currently, it
    1) produces a .seqs file of all designed sequences,
    2) optionally produces a "strands to order" file and
    3) tests the specified kinetic paths
  
  In the future it might test for bad kinetics, etc.
  """
  
  print "Finishing compilation of %s ..." % savename
  # Re-load the design system/component
  system = load(savename)
  
  print "Applying the design from '%s'" % designname
  seqs = read_design(designname)
  apply_design(system, seqs)
  
  # Document all sequences, super-sequences, strands, and structures
  print "Writing sequences file: %s" % seqsname
  f = open(seqsname, "w")
  
  f.write("# Sequences\n")
  for name, seq in system.base_seqs.items():
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
  if strandsname:
    print 'Writing "stands to order" file: %s' % strandsname
    f = open(strandsname, "w")
    for name, strand in system.strands.items():
      if not strand.dummy:
        f.write("strand %s\t%s\n" % (name, strand.seq) )
    f.close()
  
  # TODO: Write a thermodynamic scorecard.
  
  # Run Kinetic tests
  if run_kin:
    print "Testing Kinetics with parameters: ", keys
    for gate, prefix in get_gates(system):
      kinetic(gate, prefix, cleanup, **keys)


def apply_design(system, seqs):
  """Assigns designed sequences and provided mfe structures to the respective objects."""
  # Assign all the designed sequences
  for name, seq in system.base_seqs.items():
    seq.seq  = seqs[name]
    assert len(seq.seq) == seq.length
    seq.wc.seq = wc(seq.seq)
    assert seq.wc.seq == seqs[name + "*"]
  
  for sup_seq in system.sup_seqs.values():
    sup_seq.seq  = string.join([seq.seq for seq in sup_seq.base_seqs], "")
    sup_seq.wc.seq = string.join([seq.seq for seq in sup_seq.wc.base_seqs], "")
  
  for strand in system.strands.values():
    strand.seq = string.join([seq.seq for seq in strand.base_seqs], "")
  
  for name, struct in system.structs.items():
    struct.seq = string.join([strand.seq for strand in struct.strands], "+")
    assert struct.seq == seqs[name], "Design is inconsistant! %s != %s" % (struct.seq, seqs[name])


# TODO: get rid of need for gate, so get rid of recursion.
def kinetic(gate, prefix, cleanup, **keys):
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
    forward, reverse, overtime, summary = test_kinetics(kin, cleanup, **keys)
    # Process results
    num_for = len(forward)
    num_rev = len(reverse)
    num_over = len(overtime)
    num_trials = num_for + num_rev + num_over
    
    if num_over > 0:
      print
      print "WARNING: %d/%d trajectories went overtime." % (num_over, num_trials)
      print "Reported statistics may be unreliable."
      print
    
    # Rate at which collisions that will eventually go forward happen
    for_coll_rate = sum([coll_rate for (time, coll_rate) in forward]) / num_trials
    
    for_rates = [1/time for (time, coll_rate) in forward]
    for_rate = stat.mean(for_rates)
    for_rate_stddev = stat.stddev(for_rates)
    
    print "* %d/%d trajectories went forward." % (num_for, num_trials)
    if num_for > 0:
      print "  Estimated Forward Collision Rate: %f /uM/s" % (for_coll_rate / 1000000)
      print "  Estimated Forward Trajectory Rate: %f (std-dev %f) /s" % (for_rate, for_rate_stddev)
      print
    
    # Rate at which collisions that will eventually reverse happen
    rev_coll_rate = sum([coll_rate for (time, coll_rate) in reverse]) / num_trials
    
    rev_rates = [1/time for (time, coll_rate) in reverse]
    rev_rate = stat.mean(rev_rates)
    rev_rate_stddev = stat.stddev(rev_rates)
    
    print "* %d/%d trajectories went back." % (num_rev, num_trials)
    if num_rev > 0:
      print "  Estimated Reverse Collision Rate: %f /uM/s" % (rev_coll_rate / 1000000)
      print "  Estimated Reverse Trajectory Rate: %f (std-dev %f) /s" % (rev_rate, rev_rate_stddev)

  
if __name__ == "__main__":
  import re
  from optparse import OptionParser
  
  # Parse command line options.
  usage = "usage: %prog [options] BASENAME"
  parser = OptionParser(usage=usage)
  parser.set_defaults(run_kin=True, cleanup=False)
  #parser.set_defaults(verbose=True)
  #parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
  # TODO: implement quiet
  #parser.add_option("-q", "--quiet", action="store_false", dest="verbose")
  parser.add_option("--save", help="Saved state file [defaults to BASENAME.save]", metavar="FILE")
  parser.add_option("--design", help="Design file [defaults to BASENAME.mfe]", metavar="FILE")
  parser.add_option("--seqs", help="Sequences output file [defaults to BASENAME.seqs]", metavar="FILE")
  parser.add_option("--strands", help="Produce a strands-to-order file", metavar="FILE")
  parser.add_option("--no-kin", action="store_false", dest="run_kin", help="Don't run kinetics")
  
  parser.add_option("--keep-temp", action="store_false", dest="cleanup", help="Keep temporary files (.st, .wc, .eq, .sp) [Default temporarily]")
  parser.add_option("--cleanup", action="store_true", dest="cleanup", help="Remove temporary files after use.")
  #TODO: parser.add_option("--kinetic", help="Custom kinetics output file, defaults to BASENAME.kin")
  (options, args) = parser.parse_args()
  
  # Get basename
  if len(args) < 1:
   parser.error("missing required argument BASENAME")
  basename = args[0]
  # Infer the basename if a full filename is given
  p = re.match(r"(.*)\.(save|mfe)\Z", basename)
  if p:
    basename = p.group(1)
  
  # Set filename defaults
  if not options.save:
    options.save = basename + ".save"
  if not options.design:
    options.design = basename + ".mfe"
  if not options.seqs:
    options.seqs = basename + ".seqs"
  
  # Eval remaining arguments
  import quickargs
  args, keys = quickargs.get_args(args[1:])
  assert not args
  
  finish(options.save, options.design, options.seqs, options.strands, options.run_kin, options.cleanup, keys)
