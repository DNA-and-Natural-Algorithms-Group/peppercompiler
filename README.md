This is the Winfree lab DNA Circut Compiler, peppercompiler, developed by Shawn
Ligocki and currently maintained by Constantine Evans.

This git repository is the *new* version of peppercompiler, which installs as a
Python package.  It is *not ready for use yet*, and at the moment should not be
considered to have a stable API.  In particular, the package name will soon
change from pepper to peppercompiler.

The old Readme is reproduced below, but please note that most of it no longer
applies.

# Old Readme

This is the Winfree lab DNA Circuit Compiler developed by Shawn Ligocki.
It assists DNA programmers in building DNA computers by providing a
templating language for specifying generic components and interfacing
with state of the art designers and kinetic simulators to create and
test sequences.

## Setup

Run

`$ python config.py`

to set up the DNA Circuit Compiler to work with your system. It helps to
have NUPACKHOME or VIENNAHOME environment variables set before you
begin.

## Basic usage

Run

`$ python compile.py circuit`

to compile “circuit.sys” or “circuit.comp” and produce the design
specification, “circuit.pil” (use --des for .des format used by Joe's
designer).

Process “circuit.pil” with a designer. There is a wrapper for Erik's
SpuriousC in “design/”. Run

`$ python design/spurious_design.py circuit`

to produce the design, “circuit.mfe”.

Then,

`$ python finish.py circuit`

to produce a list of all sequences, “circuit.seqs” and run Multistrand
on the designed sequences.

## Advanced Usage

Advanced options of compiler.py, spurious\_design.py and finish.py can
be examined by using the --help flag (e.g. python compiler.py --help).

Extended functionality:

-   You can fix certain sequences when you begin compilation

Ex:

`$ python compiler.py --fixed=circuit.fixed circuit`

-   You can get verbose progress information while designing with
    spuriousC

Ex:

`$ python design/spurious_design.py -v circuit`

-   finish.py can output a list of “strands to order”

Ex:

`$ python finish.py --strands=circuit.strands circuit`

Furthermore, all input and output files can be specified specifically
with options. Therefore, for example, you can run two designs on the
same system at the same time. Ex:

`$ python compiler.py --output=run1.des --save=run1.save examples/Georg_System/Circuit`\
`$ python compiler.py --output=run2.des --save=run2.save examples/Georg_System/Circuit`

(These two designs will not interfere with each other.)

## Examples

The “examples/” directory has some instructive examples.

-   examples/system1/ - Really simple dummy system without
    parametrization
-   examples/system2/ - Same example, with parametrization
-   examples/Georg\_System/ - Georg's original logic circuit
    (Nature 2006)
-   examples/Lulu/ - Some mock-ups of Lulu's seesaw systems
-   examples/Elisa/ - Self-activating and self-inhibiting
    transcriptional systems

## System/Component specification

The specifications for system and component files can be found on the
DNA wiki.

<http://dna.caltech.edu/wikis/dnawiki/index.php/DNA_compiler>

## Requirements

DNA Circuit Compiler requires:

-   Python &gt;= 2.5
-   pyparsing &gt;= 1.5.1 (currently included in the directory)
-   [Multistrand](Multistrand "wikilink") (available on the cluster at
    /research/bin/Multistrand or by CVS)
-   NUPACK mfe (available on the cluster at /research/bin/mfe or
    <http://www.nupack.org/downloads>)

design/spurious\_design.py requires:

-   spuriousC (available on the dna cluster at /research/bin/spuriousC
    or by CVS)
-   Vienna RNAfold (available on the dna cluser at
    /research/bin/RNAfold or
    <http://www.tbi.univie.ac.at/~ivo/RNA/>)

