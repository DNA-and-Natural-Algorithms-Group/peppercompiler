#! /usr/bin/env python

import string

import unittest

import component_parser
from component_parser import sequence_flag, nucleotide_flag

## Helper functions
def format_signal(signal):
  seq, struct = signal
  if struct:
    return "%s(%s)" % (seq, struct)
  else:
    return seq

def format_constraint(constraint):
  if ctype == sequence_flag:
    return cvalue
  else:
    assert ctype == nucleotide_flag
    return '"%s"' % cvalue  # Nucleotide constraints are in quotes


class TestComponentParser(unittest.TestCase):
  
  ## Declare statements tests
  example_declare = (
    ("NAME", [], [], []),
    ("cOmP_f__32", [], [], []),
    ("NAME", [], [("x", "X")], [("y", "Y")]),
    ("NAME", [], [("x", None)], [("y", None)]),
    ("NAME", [], [("x", "X"), ("cow", None), ("BullDoGG_31", "bull__dogg32Strand")], [("x", None), ("y", "X")]),
    ("NAME", ["a", "toe", "BoB"], [("x", "X")], [("y", "Y")]),
  )
  
  def test00_declare_noerror(self):
    """Test simple Component Declare statement is accepted"""
    statement = "declare component NAME: x(X) -> y(Y)"
    component_parser.parse_declare_statement(statement)
  
  def test01_declare_simple(self):
    """Test simple Component Declare statement is parsed correctly"""
    statement = "declare component NAME: x(X) -> y(Y)"
    result = component_parser.parse_declare_statement(statement)
    self.assertEqual(("NAME", [], [("x", "X")], [("y", "Y")]), result)
  
  def test02_declare_examples(self):
    """Test example Component Declare statements are parsed correctly"""
    for name, params, inputs, outputs in self.example_declare:
      # Build the statement
      params_str = string.join(params, ", ")
      inputs_str = map(format_signal, inputs)
      inputs_str = string.join(inputs_str, " + ")
      outputs_str = map(format_signal, outputs)
      outputs_str = string.join(outputs_str, " + ")
      statement = "declare component %s(%s): %s -> %s" % (name, params_str, inputs_str, outputs_str)
      # Test the statement
      result = component_parser.parse_declare_statement(statement)
      self.assertEqual((name, params, inputs, outputs), result)
  
  
  ## Sequence statement tests
  example_sequence = ( 
    ("NAME", [(nucleotide_flag, "5N")], None),
    ("cOmP_f__32", [(nucleotide_flag, "5N")], None),
    ("NAME", [(nucleotide_flag, "NSWNTC")], None), #TODO: This probably won't work
    ("NAME", [(nucleotide_flag, "5N")], 5),
    ("NAME", [(sequence_flag, "seq1")], None),
    ("NAME", [(sequence_flag, "seq1*")], None), #TODO: complementary sequences
    ("NAME", [(sequence_flag, "seq1"), (nucleotide_flag, "5N"), (sequence_flag, "seq2")], None),
  )
  
  def test10_sequence_noerror(self):
    """Test simple Component Sequence statement is accepted"""
    statement = 'sequence NAME = "5N"'
    component_parser.parse_sequence_statement(statement)
  
  def test11_sequence_simple(self):
    """Test simple Component Sequence statement is parsed correctly"""
    statement = 'sequence NAME = "5N"'
    result = component_parser.parse_sequence_statement(statement)
    # We don't test constraints, because it could be parsed as 5N or N + N + N + N + N
    self.assertEqual( ("NAME", [(nucleotide_flag, "5N")], None), result)
  
  def test12_sequence_super(self):
    """Test simple Component Sequence statement is parsed correctly"""
    statement = 'sequence NAME = seq1'
    result = component_parser.parse_sequence_statement(statement)
    # We don't test constraints, because it could be parsed as 5N or N + N + N + N + N
    self.assertEqual( ("NAME", [(sequence_flag, "seq1")], None), result)
  
  def test13_sequence_examples(self):
    """Test example Component Sequence statements are parsed correctly"""
    for name, constraints, length in self.example_sequence:
      # Build the statement
      constr_str = map(format_constraint, constraints)
      constr_str = string.join(constr_str, " ")
      length_str = (": %d" % length if length != None else "")
      statement = "sequence %s = %s %s" % (name, constr_str, length_str)
      # Test the statement
      result = component_parser.parse_sequence_statement(statement)
      self.assertEqual((name, constraints, length), result)
  
  
  ## Strand statement tests
  example_strand = ( 
    ("NAME", False, [(sequence_flag, "seq")], None),
    ("cOmP_f__32", False, [(sequence_flag, "name")], None),
    ("NAME", True, [(sequence_flag, "seq")], None),
    ("NAME", False, [(sequence_flag, "seq")], 13),
    ("NAME", False, [(sequence_flag, "seq"), (sequence_flag, "OtherSequence__42")], None),
    ("NAME", False, [(sequence_flag, "seq*")], None), #TODO: complementary sequences
    ("NAME", False, [(nucleotide_flag, "5N")], None),
    ("NAME", [(nucleotide_flag, "NSWNTC")], None), #TODO: This probably won't work
  )
  
  def test20_strand_noerror(self):
    """Test simple Component Strand statement is accepted"""
    statement = 'strand NAME = seq'
    component_parser.parse_strand_statement(statement)
  
  def test21_strand_simple(self):
    """Test simple Component Strand statement is accepted"""
    statement = 'strand NAME = seq'
    result = component_parser.parse_strand_statement(statement)
    self.assertEqual( ("NAME", False, [(sequence_flag, "seq")]), result )
  
  def test22_strand_examples(self):
    """Test example Component Strand statements are parsed correctly"""
    for name, dummy, constraints, length in self.example_strand:
      # Build the statement
      dummy_str = ("[dummy]" if dummy else "")
      constr_str = map(format_constraint, constraints)
      constr_str = string.join(constr_str, " ")
      length_str = (": %d" % length if length != None else "")
      statement = "strand %s %s = %s %s" % (dummy_str, name, constr_str, length_str)
      # Test the statement
      result = component_parser.parse_strand_statement(statement)
      self.assertEqual((name, dummy, constraints, length), result)
  
  
  ## Structure statement tests
  example_structure = (
    ("NAME", ["strand"], False, ".."),
    ("cOmP_f__32", ["strand"], False, ".."),
    ("NAME", ["strand1", "St4R__and42"], False, ".."),
    ("NAME", ["strand"], True, ".."),
    ("NAME", ["strand"], False, "..(((+))..)"),
    ("NAME", ["strand"], False, "5. 3( 7( + 7) .. 3)"), #TODO: This probably won't work
    ("NAME", ["strand"], False, "U5 H3(H7(+) U2)"), #TODO: This probably won't work
  )
  
  def test30_structure_noerror(self):
    """Test simple Component Structure statement is accepted"""
    statement = 'structure NAME = strand : ..'
    component_parser.parse_structure_statement(statement)
  
  def test31_structure_simple(self):
    """Test simple Component Structure statement is accepted"""
    statement = 'structure NAME = strand : ..'
    result = component_parser.parse_structure_statement(statement)
    self.assertEqual( ("NAME", ["strand"], ".."), result )
  
  def test32_structure_examples(self):
    """Test example Component Structure statements are parsed correctly"""
    for name, strands, domain, struct in self.example_structure:
      # Build the statement
      strands_str = string.join(strands, " + ")
      domain_str = ("domain" if domain else "")
      statement = "structure %s = %s : %s %s" % (name, strands_str, domain_str, struct)
      # Test the statement
      result = component_parser.parse_structure_statement(statement)
      self.assertEqual((name, strands, domain, struct), result)
  
  
  ## Kinetic statement tests
  # TODO: kinetic constraints
  example_kinetic = (
    (["A"], ["B"]),
    ([], ["B"]),
    (["A"], []),
    (["A"], ["B", "A", "foo3Stru__ct99"]),
  )
  
  def test40_kinetic_noerror(self):
    """Test simple Component Kinetic statement is accepted"""
    statement = 'kinetic A -> B'
    component_parser.parse_kinetic_statement(statement)
  
  def test41_kinetic_simple(self):
    """Test simple Component Kinetic statement is accepted"""
    statement = 'kinetic A -> B'
    result = component_parser.parse_kinetic_statement(statement)
    self.assertEqual( (["A"], ["B"]), result )
  
  def test42_structure_examples(self):
    """Test example Component Kinetic statements are parsed correctly"""
    for inputs, outputs in self.example_kinetic:
      # Build the statement
      inputs_str  = string.join(inputs, " + ")
      outputs_str = string.join(outputs, " + ")
      statement = "kinetic %s -> %s" % (inputs_str, outputs_str)
      # Test the statement
      result = component_parser.parse_kinetic_statement(statement)
      self.assertEqual((inputs, outputs), result)

if __name__ == '__main__':
  unittest.main()

