"""Wrapper for Vienna RNAfold."""
import os
import subprocess

import RNAfold_grammar as gram
from utils import mktemp, search_file

#Globals
RNAfold = "RNAfold"

# Find dna.par file
# If users have set the nice $VIENNAHOME environment variable, it's easy.
if "VIENNAHOME" in os.environ:
  par_file = os.path.join(os.environ["VIENNAHOME"], "dna.par")
  assert os.path.isfile(par_file)

# Otherwise search the system path
else:
  path = os.environ["PATH"].split(os.path.pathsep)
  path += [".", "~"]
  par_file = search_file("dna.par", path)

assert par_file, "DNAfold_Vienna.py cannot find your dna.par file. Please set your VIENNAHOME environmant variable to it's directory and run again."
        

BREAK = "NNNNN" # Fake sequence used for strand break
def DNAfold(seq, temp):
  """Runs Vienna RNAfold on sequence 'seq' at temp 'temp' (with partition function calc if 'pf' is True)."""
  infile, infilename  = mktemp(mode="w", prefix="rna_", suffix=".in")
  outfilename = infilename[:-3] + ".out"

  # Make dummy sequence which has our fake strand break (because RNAfold can't handle multistranded folding)
  dummy_seq = seq.replace("+", BREAK)

  # Write inputs
  infile.write("%s\n" % dummy_seq)
  infile.close()
  
  # Call RNAfold
  command = "%s -T %f -P %s < %s > %s" % (RNAfold, temp, par_file, infilename, outfilename)
  subprocess.check_call(command, shell=True)
  
  # Read results
  struct, dG = gram.parseFile(outfilename)
  
  # Clean up
  os.remove(infilename)
  os.remove(outfilename)

  # Fix the structure to account for the fake strand breaking
  while BREAK in dummy_seq:
    n = dummy_seq.find(BREAK)
    dummy_seq = dummy_seq[:n] + "+" + dummy_seq[n+len(BREAK):]
    struct = struct[:n] + "+" + struct[n+len(BREAK):]
  assert seq == dummy_seq
  assert len(struct) == len(dummy_seq)
  
  return struct, dG
