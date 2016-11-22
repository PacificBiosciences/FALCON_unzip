# FALCON_unzip

FALCON-Unzip contains the modules that works with FALCON ( https://github.com/PacificBiosciences/FALCON,  https://github.com/PacificBiosciences/FALCON-integrate) for full diploid assembly (representing haplotype specific contigs as "haplotigs" as assembly output). A manuscript is deposited in BioRxiv (http://biorxiv.org/content/early/2016/06/03/056887) for our evaluation and showing the application of the algorithms and the software on multiple diploid genomes. You can find more information here: https://github.com/PacificBiosciences/FALCON/wiki/FALCON-FALCON-Unzip-%22For-Phased-Diploid-Genome-Assembly-with-Single-Molecule-Real-Time-Sequencing%22

**Note that FALCON-Unzip is under continuous development work internally. We will be keeping improving the code and algorithm for more efficient computations and advancement on the accuracy of assembly results. We are in the process refactoring the code to enhance its maintainability.** 

**Currently, the code here is in its infant stage. If you are not familiar with working from source code "as-is" and making necessary modifications for your computational environment, you might want to wait until the code is more mature, better packaged and  integrated. In the meantime, you might be able to contact PacBio representives to discuss your interests on getting diploid assembly with PacBio data in a convenient way.** 

**If you try to execute the code and you experience problems with installation and set-ups, you could submit through GitHub, however, due to limited bandwidth, we might not be able to address the questions unless they are directly related to the development path. We are investing time and resources in developing the code and ensuring its robustness.**

For end-to-end work, see `example/unzip.sh` . You do need to have run it in a working directory with Falcon output and all original Daligner overlap data.

For generating quiver consensus, you will need PacBio raw sequence data with "signal pulse" information in `bam` format. You can configure the bam-file input using a configuration sent to `fc_quvier.py`.

For generating draft haplotig assembly, you need Falcon (please install the latest code from https://github.com/PacificBiosciences/FALCON-integrate). It also use `nucmer` for aligning contigs and finding redundant contigs. `nucmer` should be installed in your execution path.

For "quivering" the draft haplotigs, you will need PacBio SMRTanalysis Software suit. More specifically, you need the following executables:
```
samtools
pbalign
makePbi
variantCaller
```

You will also need extra a python package `pysam`.

-
J.C. Nov, 2016

