This is the Winfree lab DNA Circut Compiler, peppercompiler.  It assists DNA
programmers in building DNA computers by providing a templating language for
specifying generic components and interfacing with state of the art designers
and kinetic simulators to create and test sequences.

Peppercompiler can be used directly, through scripts starting with "pepper-", or
can be used as a Python library, though the API should not be considered stable
yet.  Users of older versions should note that script names have changed.

Peppercompiler also includes, and installs, SpuriousSSM.

## Contact

The Pepper Circuit Compiler was originally developed by Shawn Ligocki, 
with later modification by Chris Berlind, Joseph Schaeffer, Constantine Evans, and Erik Winfree.


## Setup and Installation

The easiest way to install peppercompiler is via Pip:

    pip install git+https://github.com/DNA-and-Natural-Algorithms-Group/peppercompiler.git
	
Alternatively, normal python installation methods (easy_install, setup.py) may
be used.  Peppercompiler is compatible with both Python 2.7 and Python 3.

Installation requires a working C compiler so that SpuriousSSM can be compiled.
If compilation does not work for you, please let us know.

To use Nupack, NUPACKHOME should be set to your Nupack directory.

To use Multistrand, make sure that its executables are on your path.  However, the Multistrand functionality has not recently been tested.

## Basic usage

Run

    $ pepper-compile circuit

to compile “circuit.sys” or “circuit.comp” and produce the design
specification, “circuit.pil” (use --des for .des format used by the Nupack designer).

Process “circuit.pil” with a designer.  For example, to use spuriousSSM, run

    $ pepper-design-spurious circuit

to produce the design, “circuit.mfe”.

Then, run

    $ pepper-finish circuit

to produce a list of all sequences, “circuit.seqs” and run Multistrand
on the designed sequences.

## Advanced Usage

Advanced options of `pepper-compile`, `pepper-design-spurious` and `pepper-finish` can
be examined by using the --help flag (e.g. python compiler.py --help).

Extended functionality:

-   You can fix certain sequences when you begin compilation
    
    Ex:

    `$ pepper-compile --fixed=circuit.fixed circuit`

-   You can get verbose progress information while designing with
    spuriousC
    
    Ex:

    `$ pepper-design-spurious -v circuit`

-   finish.py can output a list of “strands to order”

    Ex:

    `$ pepper-finish --strands=circuit.strands circuit`

Furthermore, all input and output files can be specified specifically
with options. Therefore, for example, you can run two designs on the
same system at the same time. Ex:

`$ pepper-compile --output=run1.des --save=run1.save examples/Georg_System/Circuit`
`$ pepper-compile --output=run2.des --save=run2.save examples/Georg_System/Circuit`

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
DNA wiki: <http://dna.caltech.edu/wikis/dnawiki/index.php/DNA_compiler>

## Changes

### 0.1.3

- Problems with spuriousSSM not always being installed properly by pip have been fixed.  
- Fixed sequences that are not used in the system now result in a warning, rather than a cryptic error.
- Failures of SpuriousSSM should now be caught and result in more useful error messages.
- A typo in the spuriousSSM energetics parameters has been fixed.
- Added `dg=` option to spuriousSSM, allowing target dGs for specified parts of sequences.