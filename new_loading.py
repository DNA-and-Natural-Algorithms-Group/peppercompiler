import re

import circuit_class
import template_class

def parse_declare(line):
    # TODO: allow comments and variable spacing. Maybe use re.VERBOSE
    p = re.match(r"declare (system|component) (\w+)(?:\((.*?)\))?: (.*) ?-> ?(.*)\n", line)
    assert p, "First line of file must declare system/component.\n" + line # TODO: add line and syntax
    type_, name, param_names, ins, outs = p.groups("")
    
    # Expand out the lists
    param_names = param_names.split()
    ins  = ins.split()
    outs = outs.split()
    
    return type_, name, param_names, ins, outs

# TODO: Make into a file-like object for efficiency and preserving original file.
def preprocess(f, params):
    """Preprocesses a file. Applies parameters and removes comments and newlines."""
    doc = []
    for line in f:
        # Strip comments
        line = re.sub(r"#.*\n", r"\n", line)
        
        # Skip empty lines
        if re.match(r"\s*\Z", line):
            continue
        
        # Replace all <expression> with eval(expression) (Caution: security hole)
        def eval_brackets(s):
            return str(eval(s.group(1), params))
        line = re.sub(r"<(.*?)>", eval_brackets, line)
        
        doc.append(line)
    return doc

def load_file(filename, args):
    """Load a system or component from filename."""
    f = open(filename, "r")
    line = f.readline()
    
    # Read first line.
    type_, name, param_names, ins, outs = parse_declare(line)
    # TODO: check that name matches filename
    assert len(args) == len(param_names), "Argument length mismatch" # TODO: more info
    
    # Do parameter substitution.
    params = dict(zip(param_names, args))
    doc = preprocess(f, params)
    
    if type_ == "system":
        return load_system(doc, name, ins, outs)
    else: # type_ == "component"
        return load_component(doc, name, ins, outs)


## System stuff
def parse_import(rest):
    # TODO: allow variable spacing. Maybe use re.VERBOSE
    p = re.match(r"(\w+)(?: as (\w+))?\n", rest)
    assert p, "Import syntax incorrect.\n" + rest # TODO: add line and syntax
    path, name = p.groups()
    
    return path, name

def parse_gate(rest):
    # TODO: allow variable spacing. Maybe use re.VERBOSE
    p = re.match(r"(\w+) = (\w+)(?:\((.*?)\))?: (.*) ?-> ?(.*)\n", rest)
    assert p, "Gate syntax incorrect" # TODO: add line and syntax
    gate_name, templ_name, templ_args, ins, outs = p.groups("")
    
    # Expand out the lists
    templ_args = templ_args.split()
    ins  = ins.split()
    outs = outs.split()
    
    return gate_name, templ_name, templ_args, ins, outs

def load_system(doc, name, ins, outs):
    """Build a system object from the commands in the file."""
    system = circuit_class.Circuit(name, None, ins, outs)
    for line in doc:
        command, rest = line.split(None, 1) # Split off first word
        if command == "import":
            path, name = parse_import(rest)
            system.add_import((path, name)) # TODO: Check rest ... and fix
        elif command == "component":
            gate_name, templ_name, templ_args, ins, outs = parse_gate(rest)
            system.add_gate(gate_name, templ_name, templ_args, ins, outs)
        else:
            assert False, "Command " + command + " not defined in system." # TODO: where
    return system


## Component stuff
def parse_seq(rest):
    # TODO: allow variable spacing. Maybe use re.VERBOSE
    p = re.match(r"(\w+) = (.*) : (\d*)\n", rest)
    assert p, "Sequence syntax incorrect.\n" + `rest` # TODO: add line and syntax
    name, const, length = p.groups("")
    
    length = int(length)
    #Process constraint list
    const_temp = const.split()
    const = []
    for item in const_temp:
        p = re.match(r"(\d+|\?)?(\w)", item)
        num, code = p.groups(1) # No number means 1
        const.append((num, code))
    
    return name, const, length


def load_component(doc, name, ins, outs):
    """Build a component object from the commands in the file."""
    # TODO: fix
    ins = [(x, None) for x in ins]
    outs = [(x, None) for x in outs]
    component = template_class.Gate(name, None, ins, outs)
    for line in doc:
        command, rest = line.split(None, 1) # Split off first word
        if command == "sequence":
            name, const, length = parse_seq(rest)
            component.add_sequence(name, const, length)
        # TODO: rest of commands
        else:
            assert False, "Command " + command + " not defined in component." # TODO: where
    return component

