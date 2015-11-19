# FALCON_unzip

See `example/unzip.sh` for end-to-end work. You do need to have run it in a working diretory with Falcon output and all original Daligner overlap data.
For generating quiver consensus, you will need PacBio raw sequence data with "signal pulse" information in `bam` format. You can config the bam-file input using a configuration sent to `fc_quvier.py`.

To use this code to generate draft haplotig assembly, it need Falcon (please install the latest code from https://github.com/PacificBiosciences/FALCON-integrate). It also use `nucmer` for aligning contigs to find redundant contigs. You should have it in your execution path.

For "quivering" the draft haplotigs, you will need PacBio SMRTanalysis Software suit. More specifically, you need the following executables:
```
samtools
pbalign
makePbi
variantCaller
```

You will also need extra a python package `pysam`.

-
J.C. 2015, 1119

