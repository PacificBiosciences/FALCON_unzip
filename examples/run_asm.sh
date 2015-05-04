fc_ovlp_filter_with_phase.py --fofn las.fofn --max_diff 40 --max_cov 40 --min_cov 1 --n_core 12 --min_len 10000 --db ../../1-preads_ovl/preads.db  --rid_phase_map ./arid_to_phase > preads.p_ovl
fc_phased_ovlp_to_graph.py preads.p_ovl --min_len 7000 > fc.log
fc_graphs_to_htigs.py --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id 000000F --rid_phase_map ./arid_to_phase --fasta preads4falcon.fasta
