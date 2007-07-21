#!/usr/bin/env python

import sys

bits = int(sys.argv[1])
filename = "%dBitAdder" % bits

f = file(filename, "w")

# imports
f.write("""# %d bit adder
import HalfAdder
import Adder
import LastAdder
import Detector
""" % bits)

# inputs
for i in range(bits):
  f.write("""
in nx%d
in  x%d
in ny%d
in  y%d""" % (i, i, i, i))

def rep(val):
  if val == 0:
    return "n"
  elif val == 1:
    return " "
  else:
    raise ValueError, repr(val)

# First gates (half adder)
f.write("\n\n")
for x in range(2): # input x bit
  for y in range(2): # input y bit
    f.write("""gate G0_%d%d = HalfAdder: %sx0 + %sy0 -> %ss0 + %sc1
""" % (x, y, rep(x), rep(y), rep((x+y)%2), rep((x+y)//2)))
# Other gates (adder)
for l in range(1,bits-1): # Layer
  f.write("\n")
  for x in range(2): # input x bit
    for y in range(2): # input y bit
      for c in range(2): # input carry bit
        f.write("""gate G%d_%d%d%d = Adder: %sx%d + %sy%d + %sc%d -> %ss%d + %sc%d
""" % (l, x, y, c, rep(x), l, rep(y), l, rep(c), l, rep((x+y+c)%2), l, rep((x+y+c)//2), l+1))
# Final gates (last adder)
if bits > 1:
  l = bits-1 # Layer
  f.write("\n")
  for x in range(2): # input x bit
    for y in range(2): # input y bit
      for c in range(2): # input carry bit
        f.write("""gate G%d_%d%d%d = LastAdder: %sx%d + %sy%d + %sc%d -> %ss%d + %sc%d
""" % (l, x, y, c, rep(x), l, rep(y), l, rep(c), l, rep((x+y+c)%2), l, rep((x+y+c)//2), l+1))

# Detectors
for l in range(bits):
  f.write("""
gate DNS%d = Detector: ns%d ->
gate  DS%d = Detector:  s%d ->""" % (l, l, l, l))

f.write("""
gate DNC%d = Detector: nc%d ->
gate  DC%d = Detector:  c%d ->
""" % (bits, bits, bits, bits))

