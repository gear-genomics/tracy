# Web Applications

Tracy features a range of companion web applications hosted at [https://www.gear-genomics.com](https://www.gear-genomics.com) for browsing trace alignments, inspecting variant calls or patching reference sequences. The web apps consist of [teal](#teal) (trace browser), [sage](#sage) (trace alignment), [indigo](#indigo) (trace decomposition and variant calling), [pearl](#pearl) (patching a reference sequence using trace information), [sabre](#sabre) (multiple trace alignment viewer) and the [Wily-DNA-Editor](#wily-dna-editor) (sequence editor).

# Teal

ToDo

# Sage

ToDo

# Indigo

[Indigo](https://www.gear-genomics.com/indigo) can be used to identify single-nucleotide variants (SNVs) and short insertions or deletions (InDels) in a Sanger Chromatogram trace. The application also supports deconvolution of heterozygous mutations: SNVs cause simple double peaks but heterozygous InDels cause a shift in the trace signal and Indigo can be used to separate the two overlapping alleles. The input screen of Indigo requires a chromatogram trace file in scf, abi, ab1 or ab format. Optionally, a left and right trimming size for this trace can be specified. We recommend using [teal](https://www.gear-genomics.com/teal) for estimating such trim sizes. Indigo also requires a reference sequence as input to identify variants. This can be either a wildtype chromatogram, a small sequence in FASTA format or a large indexed reference genome. Once these input requirements have been specified the launch analysis button kicks off tracy and the results are visualized in a separate browser tab. At the top, Indigo shows the actual trace signal. Below the trace viewer is an alignment of both deconvoluted alleles with respect to the reference and an alignment of both alleles against each other. Following the alignment, Indigo lists all identified variants including their rs identifier if it is a known polymorphism, a calling quality, estimated genotype and the basecalling position in the original trace. Please note that all variants are connected via hyperlinks to the original trace for easier browsing. At the very bottom, Indigo shows the decomposition plot. In case of heterozygous InDels you should observe two distinct minima in this plot. For instance, the provided example trace file contains a heterozygous 7bp deletion and thus, Indigo shows a minima for 0 and -7bp in the decomposition plot. All plots of Indigo can be saved in png format, zoomed and panned using the [plotly](http://help.plot.ly) buttons at the top of each chart.

# Pearl

ToDo

# Sabre

[Sabre](https://www.gear-genomics.com/sabre/) is a multiple sequence alignment
viewer which allows you to view all trace sequences in parallel. An alignment in
FASTA format has to provided via file upload or pasting into the provded text area.
The default number of characters per alignment row is set to 80 by default,
but this can be adjusted.
Once everything is ready, press the `Launch Analysis` button.
Note that as with all GEAR apps, we provide an example to view. For this, first
press `Load Example` and then `Launch Analysis`.

![Sabre input page](./img/sabre-input.png)

You will now get redirected to the output page where you can browse the multiple
sequence alignment of all traces. If you hover over individual bases, information
about this particular position (both in the context of the corresponding sequence
as well as the alignment) is shown on the right.

![Sabre output page](./img/sabre-output.png)

# Wily-DNA-Editor

ToDo
