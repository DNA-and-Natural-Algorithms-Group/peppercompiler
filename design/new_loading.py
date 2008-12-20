import re

import nupack_in_class

def load_file(filename):
    """Load a NUPACK style input file."""
    f = open(filename, "r")
    
    # Create new specification object ...
    spec = nupack_in_class.Spec()
    # ... and populate it with the contents of the file.
    for line in f:
        # Strip comments and whitespace
        line = re.sub(r"#.*\n", "", line)
        line = line.strip()
        # Skip empty lines
        if not line:
            continue
        
        #print line,
        # Read the command name off (if there is a command name)
        pieces = line.split() # TODO: Put in try/except for when this fails
        if pieces[0] == "structure":
            name, struct = parse_struct(line)
            spec.add_structure(name, struct)
        elif pieces[0] == "sequence":
            name, constraints = parse_seq(line)
            spec.add_sequence(name, constraints)
        elif pieces[0] == "prevent":
            pass
        elif pieces[1] == ":": # Applying sequences to a structure
            struct_name, seqs = parse_apply(line)
            spec.add_apply(struct_name, seqs)
        elif pieces[1] == "<": # Set the objective function for a structure
            pass
        else:
            assert False, "Command not valid.\n" + line # TODO: where
    
    return spec
            

def parse_struct(line):
    command, name, eq, struct = line.split(None, 3) # TODO: Put in try/except for when this fails
    assert eq == "=", "Structure syntax incorrect" # TODO: add line and syntax
    return name, struct

def parse_seq(line):
    sequence, name, eq, const = line.split(None, 3) # TODO: Put in try/except for when this fails
    assert eq == "=", "Sequence syntax incorrect" # TODO: add line and syntax
    
    #Process constraint list
    const_temp = const.split()
    const = []
    for item in const_temp:
        p = re.match(r"(\d+)?(\w)", item) # TODO: Catch parse error
        num, code = p.groups(1) # No number means 1
        const.append((int(num), code))
    
    return name, const

def parse_apply(line):
    struct_name, colon, seqs = line.split(None, 2)
    
    temp = seqs.split()
    seqs = []
    for item in temp:
        p = re.match(r"([\w_-]+)(\*)?", item)
        seq_name, wc = p.groups("")     # wc = "*" if this is the compliment, or "" otherwise
        seqs.append((seq_name, wc))
    
    return struct_name, seqs

if __name__ == "__main__":
    import sys
    
    load_file(sys.argv[1])

