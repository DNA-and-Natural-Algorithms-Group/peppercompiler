"""Selector for DNAfold."""
import sys

from .utils import error

import os
if 'NUPACKHOME' in os.environ:
  nupack_path = os.path.join(os.environ['NUPACKHOME'],'bin','mfe')
  dnafold_choice = 'nupack'
else:
  dnafold_choice = 'none'
    

def DNAfold(seq, temp=25):
  """Run the installed thermodynamic mfe package"""
  if dnafold_choice == 'nupack':
    from . import DNAfold_Nupack
    return DNAfold_Nupack.DNAfold(seq, temp, exe=nupack_path)
  else:
    raise Exception("Sorry Dave, I can't do that: running DNAfold requires that NUPACKHOME is set.")
