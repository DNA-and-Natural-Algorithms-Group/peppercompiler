"""Selector for DNAfold."""
import os
from subprocess import CalledProcessError

if not os.path.isfile("config_choices.py"):
  print error("DNA Circuit Compiler is not configured, please run config.py")
  sys.exit(1)

import config_choices as config

def DNAfold(seq, temp=25):
  """Run the installed thermodynamic mfe package"""
  if config.thermo == config.NUPACK:
    import DNAfold_Nupack
    print "Nupack", config.nupack_path
    return DNAfold_Nupack.DNAfold(seq, temp, exe=config.nupack_path)
  else:
    import DNAfold_Vienna
    print "Vienna", config.vienna_path, config.par_file
    return DNAfold_Vienna.DNAfold(seq, temp, exe=config.vienna_path, par_file=config.par_file)
