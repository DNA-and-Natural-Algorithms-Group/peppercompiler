"""Useful classes"""
import string

class ordered_dict(dict):
  """A standard dictionary that remembers the order you added items in.
     don't use methods clear, copy, iter*, pop*, update, fromkeys."""
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
  def _dummy(*args, **kw):
    raise Exception, "methods not available"
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

class PrintObject(object):
  """Generic default-printable object."""
  def __str__(self):
    attribs = ["%s=%r" % (name, value) for (name, value) in self.__dict__.items()]
    attribs = string.join(attribs, ", ")
    return "%s(%s)" % (self.__class__.__name__, attribs)
  __repr__ = __str__

