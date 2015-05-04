~/bin/blasr_yli_0317 000000F_reads.fa 000000F_ref.fa  -noSplitSubreads -clipping subread  -placeRepeatsRandomly -minPctIdentity 70.0  -minMatch 12  -nproc 48 -bam -out aln.bam
#~/bin/blasr_yli_0317 ../data/chm1_bam/all_merged.bam pa_ctg.fa  -noSplitSubreads -clipping subread  -placeRepeatsRandomly -minPctIdentity 70.0 -minMatch 12 -nproc 48 -bam -out test_chm1.bam 
#samtools sort -m 120G test_chm13.bam chm13_aln 
#samtools sort -m 120G test_chm1.bam chm1_aln 
#python screen_aln_chm1.py
#python screen_aln_chm13.py
#bamtools merge -in chm1_aln_sampled.bam -in chm13_aln_sampled.bam -out chm1_chm13_sampled_merged.bam
#samtools sort -m 120G chm1_chm13_sampled_merged.bam chm1_chm13_sampled_merged_sorted
#samtools index chm1_chm13_sampled_merged_sorted.bam
#variantCaller.py --skipUnrecognizedContigs -v -j40 -r pa_ctg.fa  --algorithm=quiver chm1_chm13_sampled_merged_sorted.bam -o cns.fastq.gz -o cns.fasta.gz
