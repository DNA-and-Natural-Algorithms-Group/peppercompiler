"""Designs a sequence randomly."""

import sys, random, string, re
from nupack_in_parser import load_design

def random_choice(group):
  group = tuple(group)
  assert len(group) > 0
  return random.choice(group)

# Global DNA nt groups
group = {"A": "A", "T": "T", "U": "T", "C": "C", "G": "G",
         "W": "AT", "S": "CG", "N": "ATCG"} #... Others can be put later if needed ...
compliment = {"A": "T", "T": "A", "C": "G", "G": "C"}

class default_list(list):
  def __init__(self, default_func):
    self.default = default_func
  def __getitem__(self, index):
    if len(self) <= index:
      self.extend([self.default(a) for a in range(len(self), index+1)])
    return list.__getitem__(self, index)

class cool_set(set):
  def __init__(self, *args):
    self.base = None
    set.__init__(self, *args)

def set_str(string, i, val):
  assert len(val) == 1
  return string[:i] + val + string[i+1:]

class Connect(object):
  def __init__(self):
    self.data = default_list(lambda a: default_list(lambda b: [cool_set(), cool_set([(a,b)])]))
  def add(self, compl, x, y):
    # sequence number, location on sequence, parity (answers: is not wc compliment?)
    seq_x, loc_x, par_x = compl.seq_loc(x)
    seq_y, loc_y, par_y = compl.seq_loc(y)
    ### TODO: if seq_x, loc_x, par_x == seq_y, loc_y, par_y: return  # Speedup
    data_x = self.data[seq_x][loc_x]
    data_y = self.data[seq_y][loc_y]
    # Equate all involved appropriately
    data_x[par_x].update(data_y[not par_y])  
    data_x[not par_x].update(data_y[par_y])
    for i,j in data_x[True]:
      self.data[i][j][True]  = data_x[True]
      self.data[i][j][False] = data_x[False]
    for i,j in data_x[False]:
      self.data[i][j][True]  = data_x[False]
      self.data[i][j][False] = data_x[True]
    """
    foo = (1,0)
    if foo in data_x[True]  or foo in data_x[False]:
      print compl.name, x, y
      print d.seqs.get_index(seq_x).name, d.seqs.get_index(seq_y).name
      for seq in compl.seqs:
        print seq.name,
      print
      print data_x
      raw_input()"""

def design(infilename, outfile):
  global d # DEBUG
  d = load_design(infilename)
  
  # Create complementarity matrix
  connect = Connect()
  for compl in d.complexes.values():
    for x,y in compl.bonds:
      connect.add(compl, x, y)
  # Randomly color sequences
  for i, seq in enumerate(d.seqs.values()):
    for j, symb in enumerate(seq.seq):
      data = connect.data[i][j]
      if not data[True].base:
        grp = set(group[symb])
        for (x,y) in data[True]:
          grp.intersection_update(group[d.seqs.get_index(x).seq[y]])
        for (x,y) in data[False]:
          grp.intersection_update( [compliment[symb] for symb in group[d.seqs.get_index(x).seq[y]]] )
        #print data
        #print grp
        #while True:  print input()
        data[True].base = random_choice(grp)
        data[False].base = compliment[ data[True].base ]
      
      seq.seq = set_str(seq.seq, j, data[True].base)
  
  output_sequences(d, connect, outfile)

def output_sequences(d, connect, fn):
  f = file(fn, "w")
  for seq in d.seqs.values():
    # Write sequence (with dummy content)
    f.write("%d:%s\n" % (0, seq.name))
    gc_content = 0 # (seq.seq.count("C") + seq.seq.count("G")) / seq.length
    f.write("%s %f %f %d\n" % (seq.seq, 0, gc_content, 0))
    f.write("..\n..\n")
    # Write wc of sequence (with dummy content)
    f.write("%d:%s*\n" % (0, seq.name))
    wc_seq = string.join([compliment[symb] for symb in seq.seq[::-1]], "")
    f.write("%s %f %f %d\n" % (wc_seq, 0, gc_content, 0))
    f.write("..\n..\n")
  f.write("Total n(s*) = %f" % 0)
  f.close()

in_name = sys.argv[1]
out_name = re.sub(r"\.des\Z", "", in_name) + ".summary"
design(in_name, out_name)

