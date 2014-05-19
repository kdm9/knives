knives
======


A pipelined quality control script that utilises the UC Davis Biocore's
[sabre](https://github.com/najoshi/sabre), [sickle](https://github.com/najoshi/sickle),
[scythe](https://github.com/vsbuffalo/scythe) and [seqqs](https://github.com/vsbuffalo/seqqs).

Usage:
------

Knives takes paired input FASTQ files, and outputs paired good quality reads,
and single reads whose pairs failed quality control. It also takes a range of
mandatory and optional parameters. Run

    $ knives

to see a full range of parameters. A detailed user manual is on it's way, but
as this tool is simply a wrapper around [sabre](https://github.com/najoshi/sabre),
[sickle](https://github.com/najoshi/sickle), [scythe](https://github.com/vsbuffalo/scythe)
and [seqqs](https://github.com/vsbuffalo/seqqs), one can for the time being consult
their excellent user guides, available from each tool's GitHub page.
