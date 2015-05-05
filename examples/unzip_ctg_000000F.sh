fc_track_reads.py
mkdir -p 3-falcon_unzip/000000F
cd 3-falcon_unzip/000000F
fc_fetch_reads.py --fofn ../../input.fofn --ctg_fa ../../2-asm-falcon/p_ctg.fa --ctg_id 000000F --ctg_map ../../2-asm-falcon/contig_to_read_map

blasr 000000F_reads.fa 000000F_ref.fa  -noSplitSubreads -clipping subread  -placeRepeatsRandomly -minPctIdentity 70.0  -minMatch 12  -nproc 24 -bam -out aln.bam
samtools sort aln.bam aln_sorted
samtools index aln_sorted.bam
rm aln.bam

fc_phasing.py --bam aln_sorted.bam --fasta 000000F_ref.fa --ctg_id 000000F --base_dir ../
fc_phasing_readmap.py --ctg_id 000000F --ctg_map ../../2-asm-falcon/contig_to_read_map --phased_reads phased_reads

fc_ovlp_filter_with_phase.py --fofn ../../2-asm-falcon/las.fofn --max_diff 40 --max_cov 40 --min_cov 1 --n_core 12 --min_len 10000 --db ../../1-preads_ovl/preads.db  --rid_phase_map ./rid_to_phase > preads.p_ovl
fc_phased_ovlp_to_graph.py preads.p_ovl --min_len 7000 > fc.log
fc_graphs_to_htigs.py --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id 000000F --rid_phase_map ./rid_to_phase --fasta ../../1-preads_ovl/preads4falcon.fasta
