"""Multistrand wrapper"""
from __future__ import division

import string, os, subprocess, math, time, random

STOP_FLAG = "Stop_Flag"
def DNAkinfold(strands, start_struct, stop_struct, trials, sim_time, temp, conc, num_proc=1, out_interval=-1):
  """strands = dict of strand_name : sequence used in structs
     *_struct = list of complexes (each of which is a list of strand_names and a 2ndary struct)
     temp = temperature (deg C)   conc = concentration
     sim_time = max time of sim   trials = number of trials
     num_proc = number of processes to start (trials devided between them)."""
  # Assumes that structures are connected complexes
  assert start_struct != stop_struct
  trials_each = int(math.ceil(trials / num_proc))
  trials = trials_each * num_proc
  # TODO: use a standard method to make non-coliding names in /tmp once Multistrand accepts paths.
  t = time.time()
  r = random.random()
  in_name  = ".tmp.multi.in_%r-%r" % (t, r)
  out_name = ".tmp.multi.out_%r-%r" % (t, r)
  
  # Print input file for Multistrand
  f = file(in_name, "w")
  # Strand Definitions
  f.write("#Strands\n")
  for name, seq in strands.items():
    f.write("strand_%s,%s\n" % (name, seq))
  # Start Structure
  f.write("#StartStructure\n")
  for compl in start_struct:
    struct = compl.struct.replace("+", "_")
    struct = struct.replace("_", ".") # TODO: leave "_" as soon as that works
    s = ["strand_"+x for x in compl.strands]
    f.write(string.join(s, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    f.write(struct + "\n")                        # ((...._..((_))..)..(_..)...)
  # Stop Structure
  f.write("#StopStructures\n")
  for compl in stop_struct:
    struct = compl.struct.replace("+", "_")
    struct = struct.replace("_", ".") # TODO: leave "_" as soon as that works
    s = ["strand_"+x for x in compl.strands]
    f.write(string.join(s, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    f.write(struct + "\n")
  f.write("TAG: %s\n" % STOP_FLAG)
  # Other params
  f.write("##Temperature=%f\n" % temp) # Currently not working
  f.write("#Concentration=%f\n" % conc)
  f.write("#SimTime=%d\n" % sim_time)
  f.write("#NumSims=%d\n" % trials_each)
  f.write("#Logfile=%s\n" % out_name)
  f.write("#OutputInterval=%d\n" % out_interval)   # -1 = Suppress output
  f.write("#StopOption=2\n")        # Stop on stop structures defined above
  # Done
  f.close()
  
  # Run Multistrand!
  try:
    pass #os.remove(out_name)
  except OSError:
    pass # If out_name doesn't exist, we're done
  if out_interval == -1:
    command = "Multistrand > /dev/null < %s" % in_name
  else:
    command = "Multistrand < %s" % in_name
  #print command

  # Start 'num_proc' processes
  procs = []
  for i in range(num_proc):
    procs.append( subprocess.Popen(command, shell=True) )
    time.sleep(1.0)
  # Wait for them to finish
  for i in range(num_proc):
    return_code = procs[i].wait()
    if return_code != 0:
      print
      raise OSError, "Multistrand failed with status (%d)" % return_code
  #os.remove(in_name)
  
  # Read back results
  f = file(out_name, "r")
  res = f.read()
  f.close()
  
  frac = res.count(STOP_FLAG) / trials
  times = []
  for line in res.split("\n"):
    if STOP_FLAG in line:
      start = line.find(STOP_FLAG) + len(STOP_FLAG)
      times.append(float(line[start:]))
  
  return frac, times, res

