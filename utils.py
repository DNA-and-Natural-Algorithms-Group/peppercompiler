"""Useful classes and functions."""
import tempfile
import os
import string

def mktemp(mode, *args, **keys):
  """Creates a temporary file. Returns the file and filename.
     Like tempfile.mkstemp except that it actually creates a file object."""
  fd, filename = tempfile.mkstemp(*args, **keys)
  file_ = os.fdopen(fd, mode)
  return file_, filename

