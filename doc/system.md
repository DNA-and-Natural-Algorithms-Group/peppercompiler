A system file (name.sys) contains a specification for the connectivity
of [components](DNA_compiler/component_syntax "wikilink") in a DNA
system. It specifies each component and and ties their signal sequences
together. A system may be used as a component in a larger system.

#### Comments

The pound sign (\#) denotes that the rest of the line is a comment
(python/sh style comments):

`# This is a comment`

Comments may appear on their own lines or after a command, like so:

`import And22  # Half adder is used in the first layer only`

#### System Declaration {#system_declaration}

In the spirit of Matlab\'s first line declarations, each system needs to
have a declaration statement at the first line of code. Syntax:

`declare system ``<name>`{=html}`: ``<inputs>`{=html}` -> ``<outputs>`{=html}

For example:

`declare system HalfAdder: x0 + y0  ->  s0 + c1`

The `<inputs>`{=html} and `<outputs>`{=html} are \'+\' separated lists
of sequences that will be constrained to components below.

These are used if this is a system being used as a component. However,
for top-level systems, we will not have inputs and outputs and so

`<name>`{=html} should be the same as the base name of the component
file (i.e. this file should be HalfAdder.comp).

#### Imports

In order to use a component file you must import it first
(python-style). Syntax:

`import ``<component name>`{=html}

For example:

`import And22`

You may import from a different directory/name than the name you use in
the file. For example:

`import Georg0711/Parallel_And22 as And22`

which would import from the component file \'Parallel\_And22.comp\' from
the directory \'Georg0711/\', but still use the name And22 in the
specification.

#### Components

This is the meat of the circuit file. You must make one statement for
each component in the DNA system specifying the input and output
sequences. The syntax is:

`component ``<name>`{=html}` = ``<component name>`{=html}`: ``<list of input sequences>`{=html}` -> ``<list of output sequences>`{=html}

for example

`component G0_01 = And22: nx0 + y0  ->  s0 + nc1`

Note: Components **may** have 0 inputs or 0 outputs (or even both, if
you can find that useful). Example:

`component DNS0 = Detector: ns0 ->`

Is a \'ns0\' detector. It will activate if ns0 is input but doesn\'t
produce any outputs for downstream gates.

#### Examples

-   [DNA compiler/Half Adder](DNA_compiler/Half_Adder "wikilink")
-   [DNA compiler/Two-Bit Adder](DNA_compiler/Two-Bit_Adder "wikilink")
