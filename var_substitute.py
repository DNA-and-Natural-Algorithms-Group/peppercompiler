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

import string
import re

def process(infilename, params):
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
  
  import circuit_parser
  import template_parser
  
  def substitute(filename, args):
    # Parse for function declaration
    if re.match(r".*\.sys\Z", filename):
      param_names = circuit_parser.decl_stat.parseFile(filename)[2]
    elif re.match(r".*\.comp\Z", filename):
      param_names = template_parser.decl_stat.parseFile(filename)[2]
    else:
      raise ValueError, "File %s is neither system nor component type." % filename
        
    params = {}
    assert len(param_names) == len(args), (param_names, args)
    for name, val in zip(param_names, args):
      params[name] = val
    return process(filename, params)

  filename = sys.argv[1]
  args = map(eval, sys.argv[2:])
  print substitute(filename, args)

