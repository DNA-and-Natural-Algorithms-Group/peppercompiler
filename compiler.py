import sys, pickle
from circuit_parser import load_circuit
from kinetics import read_nupack, test_kinetics

def compiler(infilename):
  # Read in circuit design
  circuit = load_circuit(infilename)
  # Prepare it for Zadeh's Design
  ### TODO-maybe: deal with false input issue, where output c1 and the 4 inputs c1 are all treated independent
  circuit.output_nupack(infilename+".des")
  # Call Zadeh's Design software
  save(circuit, infilename+".save1")
  ### TODO: Work out automatic way of doing this with Joe Zadeh.
  raw_input("Run %s.des in NUPACK, save the result in %s.summary and press enter to continue." % (infilename, infilename))
  # Read results
  seqs, mfe_structs = read_nupack(infilename+".summary")
  # Prepare for Schaffer's Multistrand
  for gate_name, gate in circuit.gates.items():
    for kin in gate.kinetics.values():
      # Call Multistrand instances
      ### TODO: deal with "muliple inputs" where c2 could be any of 4 strands
      res = test_kinetics(gate_name, kin, seqs, mfe_structs)
      # TODO: process results
  # TODO: Inform user of results

def save(obj, filename):
  """Save an object in case of crash or rerun."""
  f = file(filename, "w")
  pickle.dump(obj, f)
  f.close()

if __name__ == "__main__":
  compiler(sys.argv[1])

