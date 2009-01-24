import math

from DNA_classes import complement

def read(fname):
  f = open(fname, "rb")
  data = f.read()
  f.close()
  return data

def read_files(basename):
  eq = read(basename+".eq")
  eq = map(eval, eq.split())
  
  wc = read(basename+".wc")
  wc = map(eval, wc.split())
  
  st = read(basename+".st")
  st = list(st)
  
  return eq, wc, st

def test(basename):
  """Test eq, wc and st files for consistency with spuriousC"""
  eq, wc, st = read_files(basename)
  eq = [x-1 for x in eq]
  wc = [x-1 for x in wc]
  
  free_bases = 0
  
  assert len(eq) == len(wc) == len(st), "Length mismatch. ! %d == %d == %d" % (len(eq), len(wc), len(st))
  for i, (e, w, s) in enumerate(zip(eq, wc, st)):
    if s == " ":
      assert e == -1 and w == -2
    else:
      if e == i:
        if w == -2:
          free_bases += 1
        else:
          free_bases += .5
      assert e != -1 and e <= i, "Base must at least be equal to itself. eq[%d] = %d" % (i, e)
      assert st[e] == s, "Constraint mismatch. eq[%d] = %d, but st[%d] = %d != %d = st[%d]" % (i, e, i, s, st[e], e)
      assert eq[e] == eq[eq[e]] == wc[wc[e]]
      if w != -2:
        assert st[w] == complement[s], "Constraint mismatch. wc[%d] = %d, but ~st[%d] = %d != %d = st[%d]" % (i, w, i, complement[s], st[w], w)
  
  assert free_bases == math.floor(free_bases)
  
  print "Good file. %d bases, %d free bases" % (len(eq), free_bases)
  return True

if __name__ == "__main__":
  import sys
  import re
  
  basename = sys.argv[1]
  test(basename)
