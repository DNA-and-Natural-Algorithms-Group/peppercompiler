import math
import copy

def mean(xs):
  """Sample mean."""
  if len(xs) > 0:
    return sum(xs) / len(xs)
  else:
    return None

def var(xs):
  """Sample variance."""
  if len(xs) > 1:
    m = mean(xs)
    xs2 = [(x-m)**2 for x in xs]
    return sum(xs2) / (len(xs) - 1)
  else:
    return None

def stddev(xs):
  """Sample standard deviation."""
  v = var(xs)
  if v != None:
    return math.sqrt(v)
  else:
    return None


def median(xs):
  xs = copy.copy(xs)
  xs.sort()
  
  mid = int(len(xs) / 2)
  if len(xs) % 2 == 1: # len is odd
    return xs[mid]
  else: # len is even
    return (xs[mid-1] + xs[mid]) / 2
