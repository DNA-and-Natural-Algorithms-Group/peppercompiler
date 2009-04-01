"""Multistrand wrapper"""
from __future__ import division

import string
import os
import subprocess

from utils import mktemp, urandom

def random_seed():
  """
  Return a random seed x uniform random [0 <= x < 2**31].
  We use urandom to get a more random result than random.getrandbits.
  """
  return hex(urandom.getrandbits(31)).rstrip("L")

# Constant for use in Multistrand
STOP_FLAG = "Stop_Flag"

def DNAkinfold(strands, start_struct, back_struct, stop_struct, trials_each, sim_time, temp, conc, num_proc=1, out_interval=-1):
  """
  strands = dict of strand_name : sequence used in structs
  *_struct = list of structures (each of which is a list of strand_names and a secondary struct)
  temp = temperature (deg C)
  conc = concentration (mM)
  sim_time = max time of sim (approx seconds)
  num_proc = number of processes to start
  trials_each = number of trials per process
  """
  # Note: Assumes that structures are connected complexes
  assert start_struct != stop_struct
  
  trials = trials_each * num_proc
  f, in_name = mktemp(mode="w", prefix="multi_", suffix=".in")
  out_name = in_name[:-3] + ".out"
  
  # Print input file for Multistrand
  # Strand Definitions
  f.write("#Strands\n")
  for name, seq in strands.items():
    f.write("strand_%s,%s\n" % (name, seq))
  # Start Structure
  f.write("#StartStructure\n")
  for compl in start_struct:
    s = ["strand_"+x for x in compl.strands]
    f.write(string.join(s, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    struct = compl.struct.replace("+", "_")
    f.write(struct + "\n")                        # ((...._..((_))..)..(_..)...)
  # Stop Structure
  f.write("#StopStructures\n")
  for compl in stop_struct:
    s = ["strand_"+x for x in compl.strands]
    f.write(string.join(s, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    struct = compl.struct.replace("+", "_")
    f.write(struct + "\n")
  f.write("TAG: %s\n" % STOP_FLAG)
  # Other params
  f.write("##Temperature=%f\n" % temp) # Currently not working
  f.write("#Concentration=%f\n" % conc)
  f.write("#SimTime=%d\n" % sim_time)
  f.write("#NumSims=%d\n" % trials_each)
  f.write("#Logfile=%s\n" % out_name)
  f.write("#OutputInterval=%d\n" % out_interval)   # -1 = Suppress output
  f.write("#StopOptions=2\n")        # Stop on stop structures defined above
  # Done
  f.close()
  
  # Run Multistrand!
  # Start 'num_proc' processes
  procs = []
  for i in range(num_proc):
    command = "echo '#Startseed=%s' | cat - %s | nice Multistrand" % (random_seed(), in_name)
    # If we asked for quiet, keep it quiet.
    if out_interval == -1:
      command += " > /dev/null"
    print command
    procs.append( subprocess.Popen(command, shell=True) )
    #time.sleep(1.0)
  # Wait for them to finish
  for proc in procs:
    return_code = proc.wait()
    if return_code != 0:
      print
      raise OSError, "Multistrand failed with status (%d)" % return_code
  
  # Read back results
  f = open(out_name, "r")
  res = f.read()
  f.close()
  
  frac = res.count(STOP_FLAG) / trials
  times = []
  for line in res.split("\n"):
    if STOP_FLAG in line:
      start = line.find(STOP_FLAG) + len(STOP_FLAG)
      times.append(float(line[start:]))
  
  # Clean up
  #os.remove(in_name)
  #os.remove(out_name)
  
  return frac, times, res
