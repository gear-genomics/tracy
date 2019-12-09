# Usage

Tracy uses subcommands for [basecalling](#basecalling-a-chromatogram-trace-file), [alignment](#trace-alignment), [deconvolution](#deconvolution-of-heterozygous-mutations), [variant calling](#variant-calling). These subcommands are explained below.

## Basecalling a Chromatogram Trace File

[Tracy](https://github.com/gear-genomics/tracy) can basecall a trace file and output the primary sequence (highest peak) in FASTA or FASTQ format.

```bash
tracy basecall -f fasta -o out.fasta input.ab1
tracy basecall -f fastq -o out.fastq input.ab1
```

You can also get the full trace information, including primary and secondary basecalls for heterozygous variants in JSON format or as a tab-delimited text file.

```bash
tracy basecall -f tsv -o out.tsv input.ab1
tracy basecall -f json -o out.json input.ab1
```

## Trace alignment

Tracy supports a profile-to-sequence dynamic programming alignment of a trace file to a small reference nucleotide sequence in FASTA format.

```bash
tracy align -o outprefix -r reference.fa input.ab1
```

Similarly, tracy can align a trace file to a wildtype chromatogram using profile-to-profile dynamic programming.

```bash
tracy align -o outprefix -r wildtype.ab1 input.ab1
```

Alignment of a trace file to a large reference genome requires a pre-built index on a bgzip compressed genome. For instance, to align a trace file
to the human reference genome version hg38 you first need to index the reference FASTA file.

```bash
tracy index -o hg38.fa.fm9 hg38.fa.gz
samtools faidx hg38.fa.gz
```

Once that index has been created you can then align a Chromatogram trace file to the indexed genome.

```bash
tracy align -r hg38.fa.gz input.ab1
```

Obviously, the index needs to be built only once.

## Genome index

Pre-built genome indices for commonly used reference genomes are available for download [here](https://gear.embl.de/data/tracy/).


## Deconvolution of heterozygous mutations

Double-peaks in the chromatogram trace can cause alignment issues. Tracy supports deconvolution of heterozygous variants into two separate alleles.
Each allele is then aligned separately to the reference sequence.

```bash
tracy decompose -r hg38.fa.gz -o outprefix input.ab1
cat outprefix.align1 outprefix.align2
```

Tracy also supports the use of a wildtype chromatogram for decomposition or a small FASTA sequence.

```bash
tracy decompose -r wildtype.ab1 -o outprefix mutated.ab1
tracy decompose -r sequence.fa -o outprefix mutated.ab1
```

## Variant Calling

Tracy can call and annotate single-nucleotide variants (SNV) and insertions & deletions (InDels) with respect to a reference genome. For instance, to call variants on a Chromatogram trace file using the hg38 human reference genome.

```bash
tracy decompose -v -a homo_sapiens -r hg38.fa.gz -o outprefix input.ab1
```

This command also annotates rs identifiers of known polymorphisms and produces a variant call file in binary BCF format.
This output file in BCF format can be converted to VCF using [bcftools](https://github.com/samtools/bcftools).

```bash
bcftools view outprefix.bcf
```

## Using forward and reverse ab1 files to improve variant calling

If you do have forward and reverse trace files for the same expected genomic variant you can merge variant files and check consistency of calls and genotypes. Forward trace decomposition:

```bash
tracy decompose -o forward -a homo_sapiens -r hg38.fa.gz forward.ab1
```

Reverse trace decomposition:

```bash
tracy decompose -o reverse -a homo_sapiens -r hg38.fa.gz reverse.ab1
```

Left-alignment of InDels:

```bash
bcftools norm -O b -o forward.norm.bcf -f hg38.fa.gz forward.bcf
bcftools norm -O b -o reverse.norm.bcf -f hg38.fa.gz reverse.bcf
```

Merging of normalized variant files:

```bash
bcftools merge --force-samples forward.norm.bcf reverse.norm.bcf
```
