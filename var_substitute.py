"""
Does some simple processing:
 * <expr>  are evaluated
 * {a,b} are seperated (as in shell scripting)

So, the string:
 I have {two,three,four} apples worth $<(3-5+17)/4>.
goes to
 I have two apples worth $3.
 I have three apples worth $3.
 I have four apples worth $3.
"""

import re

def process(params, infilename):
  if not params.has_key("__builtins__"):
    params["__builtins__"] = None
  f_in = file(infilename, "r")
  out = ""
  
  for line in f_in:
    # Execute !! lines
    #if line[:2] == "!!":
    #  exec line[2:].lstrip()
    #  continue
    
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
        
        out = ""
        for op in ops:
          out += duplicate(start+op+end)
        return out
      else:
        return line
    out += duplicate(line)
  
  f_in.close()
  return out

