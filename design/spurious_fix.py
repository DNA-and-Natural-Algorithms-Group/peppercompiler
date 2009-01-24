import string

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

def write(fname, data):
  f = open(fname, "wb")
  f.write(data)
  f.close()
  return

def write_files(basename, eq, wc, st):
  eq = string.join(map(repr, eq), " ")
  write(basename + ".eq", eq)
  
  wc = string.join(map(repr, wc), " ")
  write(basename + ".wc", wc)
  
  st = string.join(st, "")
  write(basename + ".st", st)
  return

def test(basename):
  """Test eq, wc and st files for consistency with spuriousC"""
  eq, wc, st = read_files(basename)
  eq = [x-1 for x in eq]
  wc = [x-1 for x in wc]
  
  assert len(eq) == len(wc) == len(st), "Length mismatch. ! %r == %r == %r" % (len(eq), len(wc), len(st))
  for i, (e, w, s) in enumerate(zip(eq, wc, st)):
    if e == -1:
      st[i] = " "
      wc[i] = -2
    else:
      assert s != " ", "Constraint must be right. st[%r] = %r" % (i, s)
      assert st[e] == s, "Constraint mismatch. eq[%r] = %r, but st[%r] = %r != %r = st[%r]" % (i, e, i, s, st[e], e)
      if w != -1:
        assert st[w] == complement[s], "Constraint mismatch. wc[%r] = %r, but ~st[%r] = %r != %r = st[%r]" % (i, w, i, complement[s], st[w], w)
  
  print "Fixed"
  eq = [x+1 for x in eq]
  wc = [x+1 for x in wc]
  write_files(basename, eq, wc, st)
  return True

if __name__ == "__main__":
  import sys
  import re
  
  basename = sys.argv[1]
  test(basename)

