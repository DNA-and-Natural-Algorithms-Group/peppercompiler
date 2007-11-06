"""Wrapper for NUPACK mfe."""
import os, tempfile

import nupack_mfe_grammar as gram
import mySubprocess as subprocess

#Globals
nupack_mfe = "mfe"

### TODO: This is not the most efficient?
def tempfilename(*args, **keys):
  fd, filename = tempfile.mkstemp(*args, **keys)
  os.close(fd)
  return filename

def DNAfold(seq, temp):
  """Runs NUPACK mfe on sequence 'seq' at temperature 'temp'."""
  infilename  = tempfilename(prefix="mfe_", suffix=".in")
  prefix = infilename[:-3]
  outfilename = prefix+".mfe"

  seqs = seq.split("+")

  # Write inputs
  infile = file(infilename, "w")
  infile.write("%d\n" % len(seqs))
  for x in seqs:
    infile.write("%s\n" % x)
  for n in range(len(seqs)):
    infile.write("%d " % (n+1) )
  infile.close()
  
  # Call RNAfold
  command = "%s -T %f -material dna -multi %s" % (nupack_mfe, temp, prefix)
  #print command
  subprocess.check_call(command, shell=True)
  os.remove(infilename)
  
  # Read results
  dG, struct = gram.parseFile(outfilename)
  os.remove(outfilename)
  
  return struct, dG

