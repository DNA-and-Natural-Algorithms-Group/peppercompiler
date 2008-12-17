"""
Does some simple processing:
 * <expr>  are evaluated
 * {a,b} are seperated onto multiple lines (as in shell scripting)
 
Lines begining with #! may be used to set variables

So, the string:
 I have {two,three,four} apples worth $<(3-5+17)/4>.
goes to
 I have two apples worth $3.
 I have three apples worth $3.
 I have four apples worth $3.
"""

import string, re

def process(params, infilename):
  if "__builtins__" not in params:
    params["__builtins__"] = None
  f_in = file(infilename, "r")
  out = ""
  
  for line in f_in:
    # Execute #! lines
    if line[:2] == "#!":
      m = string.split(line[2:], "=")
      val = eval(m[-1], params)
      del m[-1]
      for x in m:
        params[x.strip()] = val
      out += line
      continue
    
    def eval_brackets(s):
      return str(eval(s.group(1), params))
    # Replace all <expression> with eval(expression) (Caution: security hole)
    line = re.sub(r"<([^<>]*?)>", eval_brackets, line)
    
    def duplicate(line):
      # Replace first {...,expr,...} with expr for each expression in list
      m = re.search(r"{[^{}]*?}", line)
      if m:
        ops = m.group()[1:-1].split(",")
        start = line[:m.start()]
        end = line[m.end():]
        this_out = ""
        for op in ops:
          this_out += duplicate(start+op+end)
        return this_out
      else:
        return line
    out += duplicate(line)
  
  f_in.close()
  return out

if __name__ == "__main__":
  import sys
  print process({}, sys.argv[1])

