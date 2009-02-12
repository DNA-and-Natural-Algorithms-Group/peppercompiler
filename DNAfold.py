"""Selector for DNAfold."""

from mySubprocess import CalledProcessError
import DNAfold_Nupack, DNAfold_Vienna

def DNAfold(seq, temp=25):
  """Run the installed thermodynamic mfe package"""
  try:
    return DNAfold_Nupack.DNAfold(seq, temp)
  except CalledProcessError, e1:
    try:
      return DNAfold_Vienna.DNAfold(seq, temp)
    except CalledProcessError, e2:
      raise Exception, "NUPACK mfe and Vienna RNAfold both failed.\n" \
                       "NUPACK with status %s\n" \
                       "RNAfold with status %s" % (e1, e2)
