#!/usr/bin/env python
"""
Designs sequences using Winfree's SpuriousDesign/spuriousC.c algorithm.
Uses Joe Zadah's input and output formats for compatibility with compiler.
"""

import string
import subprocess

from new_loading import load_file
from DNA_classes import group, rev_group, complement

# Extend path to see compiler library
import sys
here = sys.path[0] # System path to this module.
sys.path += (here+"/..")

from DNAfold import DNAfold

# HACK
#group["_"] = ""
#rev_group[""] = "_"
#complement["_"] = "_"

def intersect_groups(x1, x2):
  g1 = group[x1]; g2 = group[x2]
  inter = set(g1).intersection(set(g2))
  assert len(inter) > 0, "System overconstrained. %s and %s cannot be equal." % (x1, x2)
  inter = list(inter)
  inter.sort()
  inter = string.join(inter, "")
  return rev_group[inter]

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
    # TODO: split strands in a structure.
    
    # seq_const[n][0] = list of (s, i) where spec.structs[s].seq[m] is spec.seqs[n]
    # seq_const[n][1] = list of (s, i) where it's the compliment
    seq_const = [([], []) for seq in spec.seqs]
    
    self.structs = spec.structs.values()
    self.seqs = spec.seqs.values()
    
    # HACK: Adding psuedo-structures for sequences
    self.num_structs = len(self.structs) + len(self.seqs)
    
    self.struct_length = [None] * self.num_structs
    
    # Create seq_const and find lengths of each structure.
    for s, struct in enumerate(self.structs):
      i = 0
      for seq in struct.seqs:
        # Add this struct, index to the appropriate set.
        seq_const[seq.num][seq.reversed].append((s, i))
        i += seq.length
      self.struct_length[s] = i
    # And for each sequence psuedo-structure
    for seq_num, seq in enumerate(self.seqs):
      s = len(self.structs) + seq_num
      seq_const[seq.num][seq.reversed].append((s, 0))
      self.struct_length[s] = seq.length
    
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
  """Create eq, wc and sequence template lists in spuriousC format."""
  # Load specification from input file.
  spec = load_file(in_name)
  
  # Load structures and subsequences into Connections object
  # and connect subsequences that are equal and complimentary.
  c = Connections(spec)
  
  # Connect regions that are constrained to helixes.
  for s, struct in enumerate(c.structs):
    for start, stop in struct.bonds:
      c.apply_wc((s, start), (s, stop))
  
  # Convert from (strand, index) format to multi_index format
  # Create the conversion function
  f = {NOTHING: NOTHING}
  eq = []
  wc = []
  st = []  # We start it as a list because python strings aren't mutable
  for s, struct in enumerate(c.structs):
    start = len(eq)
    for i in xrange(c.struct_length[s]):
      this = len(eq)
      if struct.struct[this-start] == "+":
        # Single spaces between strands.
        eq += [NOTHING]
        wc += [NOTHING]
        st += [" "]
      f[(s, i)] = len(eq)
      # Get the representatives for (s, i)
      eq_si, wc_si = c.table[(s, i)]
      # and convert them
      eq.append(eq_si)
      wc.append(wc_si)
      st.append(struct.seq[i])
    # Double spaces between structures
    eq += [NOTHING, NOTHING]
    wc += [NOTHING, NOTHING]
    st += [" ", " "]
  # The psuedo-structures.
  for seq_num, seq in enumerate(c.seqs):
    s = len(c.structs) + seq_num
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
    st += seq.seq
    st += [" ", " "]
    
  
  # Update using the conversion function
  eq = [f[x] for x in eq[:-2]]  # eq[:-2] to get rid of the [NOTHING, NOTHING] at the end
  wc = [f[x] for x in wc[:-2]]
  #st = string.join(st[:-2], "") # Finally, st should be a string.
  st = st[:-2]
  
  # Constrain st appropriately
  for i in xrange(len(eq)):
    #sys.stdout.write(st[i])
    # Skip strand breaks
    if eq[i] == NOTHING:
      continue
    if eq[i] < i:
      j = eq[i]
      #temp = st[j]
      st[j] = intersect_groups(st[j], st[i])
      #if st[j] == "_" and "_" not in (st[i], temp):
      #  print
      #  print i, j, st[i], temp
        
    if wc[i] != NOTHING and wc[i] < i:
      j = wc[i]
      #temp = st[j]
      st[j] = intersect_groups(st[j], complement[st[i]])
      #if st[j] == "_" and "_" not in (st[i], temp):
      #  print
      #  print i, j, st[i], temp
  # Propogate the changes
  for i in xrange(len(eq)):
    if eq[i] != NOTHING and eq[i] < i:
      j = eq[i]
      st[i] = st[j]
    if wc[i] != NOTHING:
      j = wc[i]
      st[i] = complement[st[j]]
  
  
  # Print SpuriousC style files
  return st, eq, wc, c

def print_list(foo, filename, format):
  """Prints a list to a file using space seperation format."""
  f = open(filename, "w")
  for x in foo:
    f.write(format % x)
  f.close()

def process_result(c, inname, outname):
  """Output sequences in Joe's format."""
  # Read spuriousC's output.
  f = open(inname, "r")
  nts = f.read()
  f.close()
  # Find and seperate structures. Sequences are stored on the last line.
  nts = nts.split("\n")[-2]
  #print repr(nts)
  seqs = string.split(nts, "  ")
  
  f = open(outname, "w")
  for s, struct in enumerate(c.structs):
    # Save the sequence
    struct.seq = seqs[s].replace(" ", "+")
    struct.mfe_struct, dG = DNAfold(struct.seq)
    #print repr(struct.seq)
    
    # Write structure (with dummy content)
    f.write("%d:%s\n" % (0, struct.name))
    gc_content = (struct.seq.count("C") + struct.seq.count("G")) / len(struct.seq)
    f.write("%s %f %f %d\n" % (struct.seq, 0, gc_content, 0))
    f.write("%s\n" % struct.struct)   # Target structure
    f.write("%s\n" % struct.mfe_struct)  # Dummy MFE structure
  
  # HACK: we stored the sequences as structures.
  for seq_num, seq in enumerate(c.seqs):
    s = len(c.structs) + seq_num
    seq.seq = seqs[s]
    
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

def design(basename, verbose=False, extra_par=""):
  
  print "Preparing constraints files for spuriousC."
  # Prepare the constraints
  st, eq, wc, c = prepare(basename + ".des")
  
  # Fix divergent specifications
  eq = [x+1 for x in eq]
  for i, x in enumerate(wc):
    if x != -1:
      wc[i] += 1
  
  # Print them to files
  print_list(st, basename + ".st", "%c")
  print_list(eq, basename + ".eq", "%d ")
  print_list(wc, basename + ".wc", "%d ")
  
  if "_" in st:
    print "System overconstrained."
    sys.exit(1)
  
  # Run SpuriousC
  # TODO: take care of prevents.
  if verbose:
    quiet = "quiet=ALL | tee %s.out" % basename
  else:
    quiet = "quiet=TRUE > %s.out" % basename
  
  command = "spuriousC score=automatic template=%s.st wc=%s.wc eq=%s.eq %s %s" % (basename, basename, basename, extra_par, quiet)
  print command
  subprocess.check_call(command, shell=True)
  #subprocess.call(command, shell=True)
  
  print "Processing results of spruriousC."
  # Process results
  process_result(c, basename + ".out", basename + ".mfe")
  print "Done, results saved to %s.mfe" % basename

if __name__ == "__main__":
  import sys
  import re
  
  verbose = ("-v" in sys.argv)
  if verbose:
    sys.argv.remove("-v")
  
  try:
    basename = sys.argv[1]
    p = re.match(r"(.*)\.des\Z", basename)
    if p:
      basename = p.group(1)
  except:
    print "Usage: python spurious_design.py [-v] infilename.des"
    sys.exit(1)
  
  extra_par = ""
  if len(sys.argv) > 2:
    extra_par = string.join(sys.argv[2:], " ")
  
  design(basename, verbose, extra_par)
