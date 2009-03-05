#!/usr/bin/env python
import string

from compiler import load
from kinetics import read_nupack, test_kinetics
import myStat as stat

from circuit_class import load_file, Circuit
from gate_class import Gate

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
  seqs, mfe_structs = read_nupack(basename+".mfe")
  
  # Apply designer results to our system
  apply_design_rec(system, "", seqs, mfe_structs)
  
  # Get all strands that will be used in the final experiment
  print 'Writing "stands to order" file: %s.strands' % basename
  strands = get_strands_rec(system)
  f = open(basename + ".strands", "w")
  for strand_name, strand_seq in strands:
    f.write("#Strand: %s\t%s\n" % (strand_name, strand_seq) )
  f.close()
  
  # Run Kinetic tests
  print "Testing Kinetics with parameters: ", keys
  kinetic_rec(system, "", **keys)


def apply_design(gate, prefix, seqs, mfe_structs):
  """Assigns designed sequences and provided mfe structures to the respective objects."""
  # Assign all the designed sequences
  for seq in gate.nupack_seqs.values():
    seq.seq  = seqs[prefix + seq.name]
    assert len(seq.seq) == seq.length
    seq.wc.seq = seqs[prefix + seq.wc.name]
    assert len(seq.wc.seq) == seq.wc.length
  for sup_seq in gate.sup_seqs.values():
    sup_seq.seq  = string.join([seq.seq for seq in sup_seq.nupack_seqs], "")
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

