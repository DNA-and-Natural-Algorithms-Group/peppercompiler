import subprocess

class CalledProcessError(Exception): pass

def check_call(*args, **kw):
  try:
    stat = subprocess.call(*args, **kw)
  except OSError:
    raise CalledProcessError, 0
  if stat != 0:
    raise CalledProcessError, stat
  return

