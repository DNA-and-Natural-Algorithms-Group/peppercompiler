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


## Generic Objects

def _dummy(*args, **kw):
  """A method not available function."""
  raise Exception, "methods not available"

class ordered_dict(dict):
  """A standard dictionary that remembers the order you added items in.
     Only supports __setitem__, __iter__, keys, values, items and read-only ops."""
  def __init__(self):
    self.order = []
    dict.__init__(self)
  def __setitem__(self, key, value):
    if key not in self:
      self.order.append(key)
    dict.__setitem__(self, key, value)
  def __iter__(self):
    for key in self.order:
      yield key
  def keys(self):
    return self.order
  def values(self):
    return [self[key] for key in self]
  def items(self):
    return [(key, self[key]) for key in self]
  
  def __delitem__(self, key):
    self.order.remove(key)
    dict.__delitem__(self, key)
  ### TODO-maybe: impliment clear, copy, iter*, ...
  clear = copy = iteritems = iterkeys = itervalues = pop \
        = popitem = update = fromkeys = _dummy

class default_ordered_dict(ordered_dict):
  """An ordered dictionary automatically defaults to a given value on unset key.
     With the call=True flag, the default will be assumed a function and called each time."""
  def __init__(self, default=None, call=False):
    self.default = default
    self.call = call
    ordered_dict.__init__(self)
  def __getitem__(self, key):
    if not call:
      return self.get(key, self.default)
    else:
      return self.get(key, self.default())

class ordered_set(set):
  """A standard set that remembers the order you added items in.
     Only supports add, update, __iter__ and read-only ops."""
  def __init__(self):
    self.order = []
    set.__init__(self)
  def add(self, value):
    if value not in self:
      self.order.append(value)
      set.add(self, value)
  def update(self, other):
  	for value in other:
  		self.add(value)
  def __iter__(self):
    for value in self.order:
      yield value
  
  ### TODO-maybe: impliment clear, copy, ...
  clear = copy = difference = difference_update = discard = intersection \
  			= intersection_update = pop = remove = symmetric_difference \
        = symmetric_difference_update = union = _dummy

class PrintObject(object):
  """Generic default-printable object."""
  def __str__(self):
    attribs = ["%s=%r" % (name, value) for (name, value) in self.__dict__.items()]
    attribs = string.join(attribs, ", ")
    return "%s(%s)" % (self.__class__.__name__, attribs)
  __repr__ = __str__
