Joe Zadeh\'s HU notation is a notation for DNA/RNA secondary structure.

It describes the structure by considering strings of consecutive
unpaired or \"helically paired\" nucleotides.

Algorithmic Rule {#algorithmic_rule}
----------------

Starting from a dot-paren secondary structure, read left-to-right.

-   If you encounter \# dots (unpaired nts.) group them together as U\#
-   If you encounter \# open-parentheses which all close consecutively
    (i.e. they form a helix) group them as H\#
    -   You must now fill in the structure \"inside\" the helix (i.e.
        after the opening of the helix, but before it\'s closing).
        Because your structure is not pseudo-knotted, this structure is
        entirely contained between the start and end of the helix so it
        can be found by recursing this algorithm. You put this structure
        in parentheses after the H\#. e.g. H4( U3 H7(U2) ).
    -   Place the structure after the close of the helix after the
        H\#(???) statement.

Parentheses may not be necessary if there is exactly one token in the
helix region.

Converter
---------

There is a secondary structure converter on the DNA cluster which can
convert dotParen2HU and HU2dotParen. They are available in
/home/sligocki/HU2dotParen/. They require the pyparsing python module,
which is included in that directory. Syntax is very simple, to convert
from HU to dot-paren for example, use:

`./HU2dotParen "U3 H4(U4) U5"`

Quotes are necessary because whitespaces and parentheses are interpreted
specially by most shells. The result is printed to standard out:

`...((((....)))).....`

Examples
--------

### A single stranded strand {#a_single_stranded_strand}

In dot-paren:

`........`

In HU notation:

`U8`

### A hairpin {#a_hairpin}

In dot-paren:

`(((())))`

In HU notation:

`H4()`

### A hairpin with unpaired nt internal {#a_hairpin_with_unpaired_nt_internal}

In dot-paren:

`((((....))))`

In HU notation:

`H4(U4)`

### \... with ss region before and after hairpin {#with_ss_region_before_and_after_hairpin}

In dot-paren:

`...((((....)))).....`

In HU notation:

`U3 H4(U4) U5`

### Multiple strands {#multiple_strands}

:   *Note: the `+` represents strand break in this notation.*

In dot-paren:

`((((+))))  or  ((((_))))`

In HU notation:

`H4(+)`

### Multiple strands, a bigger example {#multiple_strands_a_bigger_example}

In dot-paren:

`...((+.((..+)).))..`

In HU notation:

`U3 H2(+ U1 H2(U2 +) U1) U2`

### A Full Adder Gate {#a_full_adder_gate}

In \"condensed\" dot-paren:

`30(     + 15( 29. + 14. 15( + 45) + 15) 6.`

In HU notation (spaces inserted to roughly line-up with the above):

`H15(H15(+ H15(U29 + U14 H15(+ ))) + )   U6`

Walkthrough of last example {#walkthrough_of_last_example}
---------------------------

We\'ll walk through specifying this structure in HU notation:

`30( + 15( 29. + 14. 15( + 45) + 15) 6.`

We see open-parentheses, so we will have a helix region. But even though
there are 30 \"(\"s, we do not have a 30-helix, because they do not end
consecutively. Only the first 15( end consecutively, so we have an
**H15**

`structure = H15(???) ???`

Well, now we have to work out the structure inside the helix, that is
this part:

`15( + 15( 29. + 14. 15( + 45) +`

We, once again, see opening parentheses and this time they all close at
the same time so we have another **H15**

`structure = H15(H15(???) ???) ???`

Now inside this helix, we have:

`+ 15( 29. + 14. 15( + 30)`

The strand-break \"+\" is converted as-is

`structure = H15(H15(+ ???) ???) ???`

Then we have another 15-helix

`structure = H15(H15(+ H15(???) ???) ???) ???`

Inside this helix is:

`29. + 14. 15( + 15)`

Finally we have some unpaired nucleotides, so we\'ll put those, and the
strand break into the structure we\'re evolving

`structure = H15(H15(+ H15(U29 + U14 ???) ???) ???) ???`

Next is the final 15-helix

`structure = H15(H15(+ H15(U29 + U14 H15(???) ???) ???) ???) ???`

Inside this helix is just the strand break

`structure = H15(H15(+ H15(U29 + U14 H15(+) ???) ???) ???) ???`

Now looking in the 3rd helix, there was nothing after the \"15)\" (which
closed the 4th helix):

`29. + 14. 15( + 15) `*`nothing`*

so we can close that helix:

`structure = H15(H15(+ H15(U29 + U14 H15(+) ) ???) ???) ???`

and ditto in the second helix after the \"30)\" (which closed the third
helix)

`+ 15( 29. + 14. 15( + 30) `*`nothing`*

so we\'ll close that one too

`structure = H15(H15(+ H15(U29 + U14 H15(+) ) ) ???) ???`

But, inside the first helix there is trailing information (the strand
break):

`15( + 15( 29. + 14. 15( + 45) `***`+`***

so we add it:

`structure = H15(H15(+ H15(U29 + U14 H15(+) ) ) +) ???`

And finally outside of the outer helix, at the very end are a few
unpaired nts

`30( + 15( 29. + 14. 15( + 45) + 15) `***`6.`***

so we stick them on the end:

`structure = H15(H15(+ H15(U29 + U14 H15(+) ) ) +) U6`
