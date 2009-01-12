#!/usr/bin/env python
from __future__ import division

import sys
import pickle
import re
import time

from circuit_class import load_gate, Circuit
from kinetics import read_nupack, test_kinetics

def compiler(infilename, args):
  # Read in circuit design
  circuit = load_gate(infilename, args)
  # TODO: allow circuit to be a gate
  if not isinstance(circuit, Circuit):
    print "Warning: compiling Gates are not completely supported yet."

  # Write the Zadeh-style design file
  outfile = file(filename+".des", "w")
  outfile.write("## Specification for %s compiled at: %s\n" % (infilename, time.ctime()))
  circuit.output_nupack("", outfile)
  outfile.close()

  # Save compiler state to be reloaded when designer finishes
  save(circuit, infilename+".save")
  ### TODO: Work out automatic way of running designer with Joe Zadeh.
  #raw_input("Run %s.des in NUPACK, save the result in %s.summary and press enter to continue." % (infilename, infilename))

def finish(infilename):
  circuit = load(infilename+".save")
  # Read results
  seqs, mfe_structs = read_nupack(infilename+".mfe")
  
  # Prepare for Schaffer's Multistrand
  # TODO: allow circuit to be a gate
  assert isinstance(circuit, Circuit)
  for gate_name, gate in circuit.gates.items():
    for kin in gate.kinetics.values():
      # Call Multistrand instances
      ### TODO: deal with "muliple inputs" where c2 could be any of 4 strands
      frac, times, res = test_kinetics(gate_name, kin, seqs, mfe_structs)
      if times:
        ave = sum(times) / len(times)
      else:
        ave = 0
      # TODO: process results
      print "kin", gate_name,
      for compl in kin.inputs:
        print compl.name,
      print "->",
      for compl in kin.outputs:
        print compl.name,
      print ":", frac, ave

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
  filename = sys.argv[1]
  p = re.match(r"(.*)\.save", filename)
  if p:
    finish(p.group(1)) # Finish started process
  else:
    args = map(eval, sys.argv[2:])
    compiler(filename, args) # Start process

