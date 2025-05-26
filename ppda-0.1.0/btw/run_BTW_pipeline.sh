#make collapse sequence

control_file=$1

	cd ____codes
	
	./__run_make_collapse_seqs.sh $control_file #run code to make collapse sequences

	cp -f ../_seqs_folder/0.filelist ../_collapse_seqs_folder/0.filelist #copy 0.filelist to the directory of collapse sequences
	
#run BASEML

	# perl __make_bin_dat_1_index_lists.pl #make bin_dat_1_index_list which is an input to run BASEML
	
	./__run_BASEML.sh #run BASEML
	
	
#make the concatenated original alignment (if multiple sequence files are concatenated in BASEML analysis)

	./__run_make_concatenated_original_seq.sh #make concatenated original alignment file
	

#run iterative BTW

	./__run_iterative_BTWest.sh $control_file #run iterative BTW method
