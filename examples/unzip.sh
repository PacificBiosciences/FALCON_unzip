# identify raw reads to each contig
fc_track_reads.py
# get the reads assigned to each contig from the initial *.subreads.fasta files
mkdir -p 3-unzip/reads/
fc_fetch_reads.py
# phasing the read 
fc_unzip.py fc_unzip.cfg
# doing the haplotig assembly, this part will be intergrated into fc_unzip.py
mkdir -p 3-unzip/1-hasm/
cd 3-unzip/1-hasm/
find ../0-phasing/ -name "rid_to_phase.*" | xargs cat  > rid_to_phase.all
fc_ovlp_filter_with_phase.py --fofn ../../2-asm-falcon/las.fofn --max_diff 120 --max_cov 120 --min_cov 1 --n_core 12 --min_len 2500 --db ../../1-preads_ovl/preads.db --rid_phase_map ./rid_to_phase.all > preads.p_ovl
fc_phased_ovlp_to_graph.py preads.p_ovl --min_len 2500 > fc.log
fc_graphs_to_h_tigs.py --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id all --rid_phase_map ./rid_to_phase.all --fasta ../../1-preads_ovl/preads4falcon.fasta
# prepare for quviering the haplotig
cd ../
rm all_phased_reads all_h_ctg_ids all_h_ctg_edges all_p_ctg_edges all_p_ctg.fa all_h_ctg.fa
find 0-phasing -name "phased_reads" | sort | xargs cat >> all_phased_reads
find 1-hasm -name "h_ctg_ids.*" | sort | xargs cat >> all_h_ctg_ids
find 1-hasm -name "p_ctg_edges.*" | sort | xargs cat >> all_p_ctg_edges
find 1-hasm -name "h_ctg_edges.*" | sort | xargs cat >> all_h_ctg_edges
find 1-hasm -name "p_ctg.*.fa" | sort | xargs cat >> all_p_ctg.fa
find 1-hasm -name "h_ctg.*.fa" | sort | xargs cat >> all_h_ctg.fa
cd ../
# identify raw reads to each primary contig or haplotig
fc_track_reads_htigs.py
# decompose PacBio raw reads with pulse information in bam files into individual bam file for each haplotig or primary contig
mkdir -p 4-quiver/reads/
fc_select_reads_from_bam.py input_bam.fofn
# run blasr and quiver to generate consensus for each haplotig or primary contig
fc_quiver.py fc_unzip.cfg
