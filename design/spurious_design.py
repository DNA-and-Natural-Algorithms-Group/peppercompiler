#!/usr/bin/env python
"""
Designs sequences using Winfree's SpuriousDesign/spuriousC.c algorithm.
Uses Joe Zadah's input and output formats for compatibility with compiler.
"""

import string
import subprocess

from new_loading import load_file

# SpuriousC notation for nothing (specifically, there being no wc compliment).
NOTHING = -1

def min_(foo):
  """Return the min element (or NOTHING if there are no elements)."""
  if len(foo) == 0:
    return NOTHING
  else:
    return min(foo)

def min_rep(*terms):
  """Return the minimum representative. Where NOTHING is no representative."""
  terms = [x for x in terms if x != NOTHING]
  return min_(terms)

class Connections(object):
  """
  Keeps track of which parts of dna sequence are constrained to be equal or 
  complimentary to which other parts.
  """
  def __init__(self, spec):
    # seq_const[n][0] = list of (s, i) where spec.structs[s].seq[m] is spec.seqs[n]
    # seq_const[n][1] = list of (s, i) where it's the compliment
    seq_const = [([], []) for seq in spec.seqs]
    
    self.num_structs = len(spec.structs.values())
    self.struct_length = [None] * self.num_structs
    # Build these sets.
    for s, struct in enumerate(spec.structs.values()):
      i = 0
      for seq in struct.seqs:
        # Add this struct, index to the appropriate set.
        seq_const[seq.num][seq.reversed].append((s, i))
        i += seq.length
      self.struct_length[s] = i
    
    self.table = {NOTHING: (NOTHING, NOTHING)}
    # Expand these sets for specific nucleotides ...
    for seq_num, (eq, wc) in enumerate(seq_const):
      seq = spec.seqs.get_index(seq_num)
      for nt in xrange(seq.length):
        new_eq = [(s, i + nt) for (s, i) in eq]
        rep_eq = min_(new_eq)
        new_wc = [(s, i + (seq.length - 1 - nt)) for (s, i) in wc]
        rep_wc = min_(new_wc)
        # ... and save those into the connection table.
        for (s, i) in new_eq:
          self.table[(s, i)] = (rep_eq, rep_wc)
        for (s, i) in new_wc:
          self.table[(s, i)] = (rep_wc, rep_eq)
  
  def apply_wc(self, x, y):
    """Apply WC complimentarity condition between items."""
    eq_x, wc_x = self.table[x]
    eq_y, wc_y = self.table[y]
    # x is complimentary to y so merge their representatives.
    new_eq_x = new_wc_y = min_rep(eq_x, wc_y)
    new_eq_y = new_wc_x = min_rep(eq_y, wc_x)
    
    for z, (eq_z, wc_z) in self.table.items():
      if eq_z != NOTHING:
        # If z == x (== y*), update it to be equal
        if eq_z in (eq_x, wc_y):
          self.table[z] = new_eq_x, new_wc_x
        # If z == y (== x*), update it
        elif eq_z in (eq_y, wc_x):
          self.table[z] = new_eq_y, new_wc_y
  
  def printf(self):
    print "Connections.printf()"
    for s in xrange(self.num_structs):
      for i in xrange(self.struct_length[s]):
        (eq, wc) = self.table[(s, i)]
        print (s, i), eq, wc
    print

def prepare(in_name):
  # Load specification from input file.
  spec = load_file(in_name)
  
  # Load structures and subsequences into Connections object
  # and connect sequences that are equal and complimentary.
  c = Connections(spec)
  
  #c.printf()
  
  # Connect regions that are constrained to helixes.
  for s, struct in enumerate(spec.structs.values()):
    for start, stop in struct.bonds:
      c.apply_wc((s, start), (s, stop))
  
  #c.printf()
  
  # Update the sequence constraints. TODO
  
  # Convert from (strand, index) format to multi_index format
  # Create the conversion function
  f = {NOTHING: NOTHING}
  eq = []
  wc = []
  st = []  # We start it as a list because python strings aren't mutable
  for s, struct in enumerate(spec.structs.values()):
    for i in xrange(c.struct_length[s]):
      f[(s, i)] = len(eq)
      # Get the representatives for (s, i)
      eq_si, wc_si = c.table[(s, i)]
      # and convert them
      eq.append(eq_si)
      wc.append(wc_si)
    # Double spaces between structures
    eq += [NOTHING, NOTHING]
    wc += [NOTHING, NOTHING]
    st += struct.seq
    st += "  "
  
  # Update using the conversion function
  eq = [f[x] for x in eq[:-2]]  # eq[:-2] to get rid of the [NOTHING, NOTHING] at the end
  wc = [f[x] for x in wc[:-2]]
  st = string.join(st, "").strip() # Finally, st should be a string.
  
  # Print SpuriousC style files
  return st, eq, wc
  #print_list(eq, basename + ".eq")
  #print_list(wc, basename + ".wc")

def print_list(foo, filename, format):
  """Prints a list to a file using space seperation format."""
  f = open(filename, "w")
  for x in foo:
    f.write(format % x)
  f.close()

#def output_sequences(d, connect, fn):
def read_result(spec, c, inname, outname):
  """Output sequences in Joe's format."""
  # Read spuriousC's output.
  f = open(inname, "r")
  nts = f.read()
  f.close()
  # And seperate structures.
  seqs = string.split(nts, "  ")
  
  f = open(outname, "w")
  for s, struct in enumerate(spec.structs.values()):
    # Save the sequence
    struct.seq = seqs[s]
    
    # Write structure (with dummy content)
    f.write("%d:%s\n" % (0, struct.name))
    gc_content = (struct.seq.count("C") + struct.seq.count("G")) / (len(struct.seq) - len(struct.strands)) # - len(struct.strands) because of the +'s in the struct.seq
    f.write("%s %f %f %d\n" % (struct.seq, 0, gc_content, 0))
    f.write("%s\n" % struct.dp_struct)   # Target structure
    f.write("%s\n" % struct.mfe_struct)  # MFE structure of chosen sequences
  
  # TODO: get sequence seqs set.
  for seq in spec.seqs.values():
    # Write sequence (with dummy content)
    f.write("%d:%s\n" % (0, seq.name))
    gc_content = (seq.seq.count("C") + seq.seq.count("G")) / seq.length
    f.write("%s %f %f %d\n" % (seq.seq, 0, gc_content, 0))
    f.write(("."*seq.length+"\n")*2) # Write out dummy structures.
    # Write wc of sequence (with dummy content)
    seq = seq.wc
    f.write("%d:%s\n" % (0, seq.name))
    #wc_seq = string.join([complement[symb] for symb in seq.seq[::-1]], "")
    f.write("%s %f %f %d\n" % (seq.get_seq(), 0, gc_content, 0))
    f.write(("."*seq.length+"\n")*2) # Write out dummy structures.
  f.write("Total n(s*) = %f" % 0)
  f.close()

def design(in_name, basename, verbose=False):
  # Prepare the constraints
  st, eq, wc = prepare(in_name)
  
  # Fix divergent specifications
  eq = [x+1 for x in eq]
  for i, x in enumerate(wc):
    if x != -1:
      wc[i] += 1
  
  # Print them to files
  print_list(st, basename + ".st", "%c")
  print_list(eq, basename + ".eq", "%d ")
  print_list(wc, basename + ".wc", "%d ")
  
  # Run SpuriousC
  # TODO: take care of prevents.
  if verbose:
    quiet = "quiet=ALL"
  else:
    quiet = "quiet=TRUE"
  
  command = "spuriousC score=automatic template=%s.st wc=%s.wc eq=%s.eq %s > %s.out" % (basename, basename, basename, quiet, basename)
  print command
  subprocess.check_call(command, shell=True)
  
  # TODO: Process results
  # = read_result(basename + ".out", basename + ".mfe")

if __name__ == "__main__":
  import sys
  import re
  
  verbose = ("-v" in sys.argv)
  if verbose:
    sys.argv.remove("-v")
  
  try:
    in_name = sys.argv[1]
    assert re.search(r"\.des\Z", in_name)
    basename = re.sub(r"\.des\Z", "", in_name) # Makes *.des => *
  except:
    print "Usage: python spurious_design.py [-v] infilename.des"
    sys.exit(1)

  design(in_name, basename, verbose)

