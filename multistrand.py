"""Multistrand wrapper"""
from __future__ import division

import string
import re
import os
import subprocess

from utils import mktemp, urandom

def random_seed():
  """
  Return a random seed x uniform random [0 <= x < 2**31].
  We use urandom to get a more random result than random.getrandbits.
  """
  return hex(urandom.getrandbits(31)).rstrip("L")

# Constants for use in Multistrand
BACK_FLAG = "REVERSE"
FORWARD_FLAG = "FORWARD"

def DNAkinfold(strands, start_struct, back_struct, stop_struct, trials, sim_time, temp, conc, out_interval=-1):
  """
  strands = dict of strand_name : sequence used in structs
  *_struct = list of structures (each of which is a list of strand_names and a secondary struct)
  temp = temperature (deg C)
  conc = concentration (mM)
  sim_time = max time of sim (approx seconds)
  trials = number of trials to run
  """
  # Note: Assumes that structures are connected complexes
  assert start_struct != stop_struct
  
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
  f.write("#StopStructures\n")
  # Backward Stop Structure (Where the reaction started and might fall back to)
  for compl in back_struct:
    s = ["strand_"+x for x in compl.strands]
    f.write(string.join(s, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    struct = compl.struct.replace("+", "_")
    f.write(struct + "\n")
  f.write("TAG: %s\n" % BACK_FLAG)
  # Forward Stop Structure (Where the reaction should go)
  for compl in stop_struct:
    s = ["strand_"+x for x in compl.strands]
    f.write(string.join(s, ",") + "\n")  # Seq1,Seq2,Chicken,Seq4
    struct = compl.struct.replace("+", "_")
    f.write(struct + "\n")
  f.write("TAG: %s\n" % FORWARD_FLAG)
  # Other params
  f.write("#Energymodel=NUPACK_DNA_2_3\n")
  f.write("#Temperature=%f\n" % temp)
  f.write("#Concentration=%f\n" % conc)
  f.write("#SimTime=%d\n" % sim_time)
  f.write("#NumSims=%d\n" % trials)
  f.write("#Logfile=%s\n" % out_name)
  f.write("#OutputInterval=%d\n" % out_interval)   # -1 = Suppress output
  f.write("#StopOptions=2\n")        # Stop on stop structures defined above
  f.write("#SimulationMode=1\n")     # Start with binding two structions and then look at forward and backwards rates.
  # Done
  f.close()
  
  # This'll be a mess if we leave an old outfile there.
  if os.path.isfile(out_name):
    os.remove(out_name)
  # Run Multistrand!
  command = "nice Multistrand < %s" % in_name
  # If we asked for quiet, keep it quiet.
  if out_interval == -1:
    command += " > /dev/null"
  print "$", command
  subprocess.check_call(command, shell=True)
  
  res = read_result(out_name)
  
  # Clean up
  #os.remove(in_name)
  #os.remove(out_name)
  
  f = open(out_name, "r")
  
  # Note: This will load the entire data of the file into lists (could be memory 
  #       intensive for extremely large data sets).
  forward = []
  reverse = []
  overtime = []
  summary = ""
  for line in f:
    if line[0] == "(":
      parts = line.split()
      assert len(parts) == 5, "Error: Multistrand format has changed.\n%s" % line
      coll_rate = float(parts[2])
      flag = parts[3]
      time = float(parts[4])
      if time >= sim_time:
        overtime.append( (time, coll_rate) )
      elif flag == FORWARD_FLAG:
        forward.append( (time, coll_rate) )
      elif flag == BACK_FLAG:
        reverse.append( (time, coll_rate) )
      else:
        assert False, "Error: Unexpected stop flag '%s' in Multistrand output.\n%s" % (flag, line)
    else:
      summary += line
  f.close()
  
  return forward, reverse, overtime, summary

class Result(object): pass

def read_result(filename):
  """
  Read the end of the results to get the rate constants
  
  Sample file "tail -7":
  ...
  Simulation Complete: 10000 trajectories total.
  Estimated Collision Reaction Rate: 0.001224 (/M/s)
  Estimated mean, var Collision Reaction Rate: 0.001224, 0.000000 (/M/s)
  Forward Trajectory Rate (mean, var): 0.002500, 0.000006 (/s)
  Forward Trajectories: 389, Average Time: 7.958826e+02
  Reverse Trajectory Rate (mean,var): 7.377919, 6048.763386 (/s)
  Reverse Trajectories: 9611, Average Time: 1.164428e+02
  """
  lines = subprocess.Popen("tail -7 %s" % filename, shell=True, stdout=subprocess.PIPE).stdout.readlines()
  
  res = Result()
  try:
    res.num_trials = re.match(r"Simulation Complete: (.+) trajectories total\.\n", lines[0]).group(1)
    res.num_trials = int(res.num_trials)
    
    res.coll_rate, res.coll_var = re.match(r"Estimated mean, var Collision Reaction Rate: (.+), (.+) \(/M/s\)\n", lines[2]).group(1,2)
    res.coll_rate = float(res.coll_rate)
    res.coll_var  = float(res.coll_var)
    
    res.for_rate, res.for_var = re.match(r"Forward Trajectory Rate \(mean(, var)?\): (N/A|([\d.e+-]+))(, ([\d.e+-]+))? \(/s\)\n", lines[3]).group(3, 5)
    if res.for_rate != None:
      res.for_rate = float(res.for_rate)
    if res.for_var != None:
      res.for_var  = float(res.for_var)
    
    res.for_num, res.for_mean_time = re.match(r"Forward Trajectories: (.+), Average Time: (N/A|([\d.e+-]+))\n", lines[4]).group(1, 3)
    res.for_num = int(res.for_num)
    if res.for_mean_time != None:
      res.for_mean_time = float(res.for_mean_time)
    
    res.rev_rate, res.rev_var = re.match(r"Reverse Trajectory Rate \(mean(,var)?\): (N/A|([\d.e+-]+))(, ([\d.e+-]+))? \(/s\)\n", lines[5]).group(3, 5)
    if res.rev_rate != None:
      res.rev_rate = float(res.rev_rate)
    if res.rev_var != None:
      res.rev_var  = float(res.rev_var)
    
    res.rev_num, res.rev_mean_time = re.match(r"Reverse Trajectories: (.+), Average Time: (N/A|([\d.e+-]+))\n", lines[6]).group(1, 3)
    res.rev_num = int(res.rev_num)
    if res.rev_mean_time != None:
      res.rev_mean_time = float(res.rev_mean_time)
  
  except AttributeError:
    print
    print lines
    print
    raise
  
  return res
