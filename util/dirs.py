import os

def updir(d, n):
  """Given path d, go up n dirs from d and return that path"""
  ret_val = d
  for _ in range(n):
    ret_val = os.path.dirname(ret_val)
  return ret_val

def makedir(adir):
    if not os.path.exists(adir):
        shell("mkdir -p " + adir)
