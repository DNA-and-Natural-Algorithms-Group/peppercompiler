from __future__ import division

import sys, pickle, re
from circuit_parser import load_circuit
from kinetics import read_nupack, test_kinetics

def finish(infilename, gate_name, kin_num):
  circuit = load(infilename+".save")
  # Read results
  seqs, mfe_structs = read_nupack(infilename+".summary")
  # Prepare for Schaffer's Multistrand
  gate = circuit.gates[gate_name]
  kin = gate.kinetics[kin_num]
  # Call Multistrand instances
  frac, times, res = test_kinetics(gate_name, kin, seqs, mfe_structs)
  ave =  ( sum(times) / len(times) if times else 0 )
  # TODO: process results
  print "kin", gate_name,
  for compl in kin.inputs:
    print compl.name,
  print "->",
  for compl in kin.outputs:
    print compl.name,
  print ":", frac, ave

def load(filename):
  """Load it back."""
  f = file(filename, "r")
  obj = pickle.load(f)
  f.close()
  return obj

if __name__ == "__main__":
  filename = sys.argv[1]
  p = re.match(r"(.*)\.save", filename)
  assert p
  finish(p.group(1), sys.argv[2], int(sys.argv[3])) # Finish started process

