Copy "codon_configuration_at_extant_nodes.txt" created in the step1.
Copy "joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW.txt" created in the step1.

run make_joint_prob_anc_cod_config_list_wo_configuration_with_stop_codons_ver221215.pl
	This code filters ancestral configurations with stop codon at any of ancestral nodes.
run make_joint_prob_anc_cod_config_list_only_with_single_or_no_change_per_path_ver221215.pl
	This code filters ancestral configuration which requires >1 nucleotide change within a lineage.

run make_list_of_codon_change_with_probability_for_dNdS_calculation_ver221215.pl
run make_list_of_codon_change_with_probability_S_N_for_fix_for_dNdS_calculation_ver221215.pl
run count_Dn_Ds.pl
	After running the three codes, you can get estimated fixation counts in each lineage using only codon configurations which requires <=1 nucleotide change per lineage.

run make_list_of_N_S_at_ms_and_ye_node.pl
	This code outputs number of synonymous(S) and nonsynonymous(N) sites estimated at ms and ye nodes assuming equal mutation rates among point mutations .

The estimated fixation counts and number of sites allow to estimate synoymous/nonsynonymous substitution rate, which will be used in the "codon_path" in step3.