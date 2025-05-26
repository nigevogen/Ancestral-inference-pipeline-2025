Copy "codon_configuration_at_extant_nodes.txt" created in the step1.
Copy "joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons.txt" created in the step2.

open "combine_BTW_and_codon_path_result_HY.py" and change the parameters "synonymous_rate" and "nonsynonymous_rate" based on the estimates in the step2

run combine_BTW_and_codon_path_result_HY.py
	This code runs codon_path and infer codon changes for each codon configuration.

Rename "test_out.txt" as "list_of_BTW_ancestor_prob_and_codon_path_prob.txt".

run make_joint_prob_anc_cod_config_list_only_with_single_or_no_change_per_path_ver221215.pl
	This code filters ancestral configuration which requires >1 nucleotide change within a lineage.

run make_list_of_codon_change_with_probability.pl
	This code makes list of codon changes with the probabilities for each position in input sequences and each lineage.
	Output file is "make_list_of_codon_change_with_probability.txt".

run convert_codon_to_digit_for_list_of_codon_change_ver221215.pl
	This code convert triplet in "make_list_of_codon_change_with_probability.txt" to our "codon_digit".

run make_list_of_codon_change_with_probability_for_Dm_Ds_fixations_ver221215_codon_digit.pl
run make_list_of_codon_change_with_probability_for_Dm_Ds_polymorphisms_ver221215_codon_digit.pl
run make_list_of_codon_change_with_probability_S_N_for_fix_ver221215_codon_digit.pl
run make_list_of_codon_change_with_probability_S_N_for_poly_ver221215_codon_digit.pl
	After running these four codes, you can get lists of inferred synonymous and nonsynonymous fixations and polymorphisms for Dmel and Dsim.
