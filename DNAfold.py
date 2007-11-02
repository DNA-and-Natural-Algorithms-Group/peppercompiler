"""Selector for DNAfold."""

from subprocess import *

#CalledProcessError
import DNAfold_Nupack, DNAfold_Vienna

def DNAfold(seq, temp):
  """Run the installed thermodynamic mfe package"""
  try:
    return DNAfold_Nupack.DNAfold(seq, temp)
  except CalledProcessError, e1:
    try:
      return DNAfold_Vienna.DNAfold(seq, temp)
    except CalledProcessError, e2:
      raise Exception, "NUPACK mfe and Vienna RNAfold both failed.\n%s\n%s" % (e1, e2)

