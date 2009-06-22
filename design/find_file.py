import os

class BadFilename(Exception):
  """File not found"""
  pass

def find_file(name):
  """Find the actual filename, which is either name or name+".des" """
  if os.path.isfile(name):
    return name
  elif os.path.isfile(name + ".des"):
    return name + ".des"
  else:
    raise BadFilename, "File not found: neither %s nor %s.des exist." % (name, name)
