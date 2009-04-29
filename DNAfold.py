"""Selector for DNAfold."""
from subprocess import CalledProcessError

import DNAfold_Nupack

def DNAfold(seq, temp=25):
  """Run the installed thermodynamic mfe package"""
  try:
    return DNAfold_Nupack.DNAfold(seq, temp)
  except CalledProcessError, e1:
    try:
      import DNAfold_Vienna
      return DNAfold_Vienna.DNAfold(seq, temp)
    except CalledProcessError, e2:
      raise Exception, "NUPACK mfe and Vienna RNAfold both failed.\n" \
                       "NUPACK with status %s\n" \
                       "RNAfold with status %s" % (e1, e2)
