#-----------------------------------------------------------------------------------------
#- copy the following code to pause the process
#- CODE: 
#-     echo ""; echo "Press enter to continue"; read Wait;
#-----------------------------------------------------------------------------------------
# setting --------------------------------------------------------------------------------
replicate=10							# number of BASEML replicates (choose the highest likelihood among replicates)
CDSopt=0								# 0: nucleotide; 1: CDS

BASEML_model_pos1="model3_pos1.ctl"							# name of the control file used for BASEML (codon pos1)
BASEML_model_pos2="model3_pos2.ctl"							# name of the control file used for BASEML (codon pos2)
BASEML_model_pos3="model3_pos3.ctl"							# name of the control file used for BASEML (codon pos3)

BASEML_model_nucleotide="model3.ctl"							# name of the control file used for BASEML (nucleotide)

bin_list_folder="sample_concat_list_codons_nucleotides"							# concat list folder
seq_folder="sample_sequences"							# sequence folder
filteropt=1						# 1: allcod; 2: 2f10cDcRY; 3, 4f6cN; 4: 2f10cDcAA; 5: 4f6cAA; 101-122: singleAA
rmopt=1									# temp folder removing option
subopt=(1 1 1)							# subset choice option (nh0, nh1 or nh4)
#-----------------------------------------------------------------------------------------
cd ../../..

rm -r 10.input
rm -r 20.output
mkdir 04_output

cd 00_bin
mkdir ../04_output/sample_result
mkdir ../10.input

for((i=1;i<=1;i++))
do
	# make output folder (folder name will be changed later) #
	mkdir ../04_output/tmp_result
	#-------------------------------------------------------------------------------------
	# make concat ctl file (for Proc2b)
	#-------------------------------------------------------------------------------------
	cd ../02_ctlfiles/sample_ctl/    # cd ctl file folder
	k="_index_lists"                          # suffix of cc_list filename
	file="bin_dat_$i$k"                   # cc_list filename
	c="2b_33_mfa_concat_$i.ctl"         # ctl filename
	echo "/*--------------------------------------------------*/" >$c
	case "$CDSopt" in
		0) echo "process_choice	= 331	/* */" >>$c;;
		1) echo "process_choice	= 331	/* */" >>$c;;
	esac		
	echo "MainFolder	= \"../04_output/tmp_result\"	/* */" >>$c
	echo "CreateFldOpt	= 1	/* */" >>$c
	echo "bin_gn_dat_list_opt     = 1     /* */" >>$c
	echo "ConcatList	= \"../../03_seq_dat/sample_data/$bin_list_folder/$file\"	/* */" >>$c
	echo "MfaFld_in		= \"../../03_seq_dat/sample_data/$seq_folder\"	/* */" >>$c
	echo "MfaFld_out	= \"2_alncds_cc.MFA\"	/* */" >>$c
	echo "NumErrorCodesOK	= 6	/* */" >>$c
	echo "ErrorCodesOK	= 0,1,2,4,8,16	/* */" >>$c
	echo "rootInputFolder	: \"../10.input/\"	/* */" >>$c
	echo "rootOutputFolder	: \"../20.output\"	/* */" >>$c
	echo "dirSep 			: \"/\"				/* */" >>$c
	echo "/*--------------------------------------------------*/" >>$c
	echo "fileOpenMode		: TEXT				/* */" >>$c
	echo "newLine			: LF				/* */" >>$c
	echo "overWriteOption	: dontOverWrite		/* */" >>$c
	echo "/*--------------------------------------------------*/" >>$c
	
	#-------------------------------------------------------------------------------------
	# run subset (nh0, Proc2b ~ Proc4)
	#-------------------------------------------------------------------------------------
	cd ../../00_bin
	#---------------------------------------
	# Proc2b: mfa -> concat_mfa
	#---------------------------------------
	b="temp_b" # batchfile for each input file
	echo "../02_ctlfiles/sample_ctl/$c" >$b
	./00_mfa_convert_4_LLVMgcc_120719 < $b
	rm -r ../20.output
	rm $b
	case "$CDSopt" in
		0) mv ../04_output/tmp_result/2_alncds_cc.MFA/0.bin_gene_dat_list ../04_output/tmp_result/0.bin_gene_dat_list;;
		1) mv ../04_output/tmp_result/2_alncds_cc.MFA/0.bin_gene_dat_list ../04_output/tmp_result/0.bin_gene_dat_list;;
	esac		

	#-----------------------------------------------------------
	# Proc2c: mfa -> selected mfa
	#---------------------------------------
	b="temp_b" # batchfile for each input file
	echo "../02_ctlfiles/sample_ctl/2c_36_mfaSelect.ctl" >$b
	./00_mfa_convert_4_LLVMgcc_120719 < $b
	rm -r ../20.output
	rm $b
	#-----------------------------------------------------------
	# Proc3: mfa -> PAML
	#---------------------------------------
	b="temp_b" # batchfile for each input file
	echo "../02_ctlfiles/sample_ctl/3_11_mfa2PAML.ctl" >$b
	./00_mfa_convert_4_LLVMgcc_120719 < $b
	rm -r ../20.output
	rm $b
	# remove mfa folder #
	if [ $rmopt -eq 1 ]; then rm -r ../04_output/tmp_result/2_alncds_cc.MFA; fi;
	
	if [ $CDSopt -eq 1 ]; then		# CDS
		#-----------------------------------------------------------
		# Proc4: PAML -> 3position PAML
		#---------------------------------------
		b="temp_b" # batchfile for each input file
		echo "../02_ctlfiles/sample_ctl/4_19_PAML3pos.ctl" >$b
		./00_mfa_convert_4_LLVMgcc_120719 < $b
		rm -r ../20.output
		rm $b
	fi;

	#-------------------------------------------------------------------------------------
	# run subsets (Proc5, Proc5b)
	#-------------------------------------------------------------------------------------
	cd ../04_output
	# make PAML output folder #
	mkdir tmp_result/5_alncds_rst_bmlo
	mkdir tmp_result/5_alncds_filter_rst
	#------------------------------------------------------------------------------
	# BASEML analysis
	#---------------------------------------------------------------	
	if [ ${subopt[0]} -eq 1 ]; then
		cd ../01_batchfiles/batch_A/batch_B
		if [ $CDSopt -eq 1 ]; then		# CDS
			#-----------------------------------------------------------
			# Proc5: BASEML for specified model
			# pos1
			#---------------------------------------
			for ((j=1;j<$replicate+1;j++)) {
			 b="batch_tmp"
			 echo "../../../02_ctlfiles/sample_ctl/5_Baseml_Ctl/$BASEML_model_pos1" >$b
			 echo "" >>$b
			
		     baseml_49f < $b
		     
			 perl check_lnL.pl > X

              filename="X"
              cat $filename | while read line
              do
                if [ "$line" = "1" ]
                 then
                 echo " "
                 echo "Higher lnL value than previous"
                 #cp highestL highestL_"$j"
			     mv -f rst ../../../04_output/tmp_result/5_alncds_rst_bmlo/pos1_rst
			     mv -f mlf ../../../04_output/tmp_result/mlf_pos1_bin"$i"
                fi;
                if [ "$line" = "0" ]
                 then
                 rm rst
                 rm mlf
                 echo " "
                 echo "Lower lnL value than previous"
                fi;
             done
 
			rm X
			rm $b
			rm 2base.t
			rm in.basemlg
			rm lnf
			rm rates			
			rm rst1
			rm rub
			 
			}
			rm highestL
			#---------------------------------------
			# pos2
			#---------------------------------------
			for ((j=1;j<$replicate+1;j++)) {
			 b="batch_tmp"
			 echo "../../../02_ctlfiles/sample_ctl/5_Baseml_Ctl/$BASEML_model_pos2" >$b
			 echo "" >>$b
			
		     baseml_49f < $b
		     
			 perl check_lnL.pl > X

              filename="X"
              cat $filename | while read line
              do
                if [ "$line" = "1" ]
                 then
                 echo " "
                 echo "Higher lnL value than previous"
                 #cp highestL highestL_"$j"
			     mv -f rst ../../../04_output/tmp_result/5_alncds_rst_bmlo/pos2_rst
			     mv -f mlf ../../../04_output/tmp_result/mlf_pos2_bin"$i"
                fi;
                if [ "$line" = "0" ]
                 then
                 rm rst
                 rm mlf
                 echo " "
                 echo "Lower lnL value than previous"
                fi;
             done
 
			rm X
			rm $b
			rm 2base.t
			rm in.basemlg
			rm lnf
			rm rates			
			rm rst1
			rm rub
			 
			}
			rm highestL
			#---------------------------------------
			# pos3
			#---------------------------------------
			for ((j=1;j<$replicate+1;j++)) {
			 b="batch_tmp"
			 echo "../../../02_ctlfiles/sample_ctl/5_Baseml_Ctl/$BASEML_model_pos3" >$b
			 echo "" >>$b
			
		     baseml_49f < $b
		     
			 perl check_lnL.pl > X

              filename="X"
              cat $filename | while read line
              do
                if [ "$line" = "1" ]
                 then
                 echo " "
                 echo "Higher lnL value than previous"
                 #cp highestL highestL_"$j"
			     mv -f rst ../../../04_output/tmp_result/5_alncds_rst_bmlo/pos3_rst
			     mv -f mlf ../../../04_output/tmp_result/mlf_pos3_bin"$i"
                fi;
                if [ "$line" = "0" ]
                 then
                 rm rst
                 rm mlf
                 echo " "
                 echo "Lower lnL value than previous"
                fi;
             done

			rm X
			rm $b
			rm 2base.t
			rm in.basemlg
			rm lnf
			rm rates			
			rm rst1
			rm rub			 
			}
			rm highestL
			
			# BASEML end #
			cd ../../../00_bin
			#-----------------------------------------------------------
			# Proc5b: BASEML output .rst -> filtered rst
			# pos1
			#---------------------------------------
		b="temp_b" # batchfile for each input file
			echo "../02_ctlfiles/sample_ctl/5b_34_rstfilter_pos1.ctl" >$b
			./00_mfa_convert_4_LLVMgcc_120719 < $b
			rm -r ../20.output
			rm $b
			#---------------------------------------
			# pos2
			#---------------------------------------
			b="temp_b" # batchfile for each input file
			echo "../02_ctlfiles/sample_ctl/5b_34_rstfilter_pos2.ctl" >$b
			./00_mfa_convert_4_LLVMgcc_120719 < $b
	        rm -r ../20.output
			rm $b
			#---------------------------------------
			# pos3
			#---------------------------------------
			b="temp_b" # batchfile for each input file
			echo "../02_ctlfiles/sample_ctl/5b_34_rstfilter_pos3.ctl" >$b
			./00_mfa_convert_4_LLVMgcc_120719 < $b
			rm -r ../20.output
			rm $b
		fi;
		if [ $CDSopt -eq 0 ]; then		# nucleotide
			#-----------------------------------------------------------
			# Proc5: BASEML for specified model
			#---------------------------------------
			for ((j=1;j<$replicate+1;j++)) {
			b="batch_tmp"
			echo "../../../02_ctlfiles/sample_ctl/5_Baseml_Ctl/$BASEML_model_nucleotide" >$b
			echo "" >>$b
			baseml_49f < $b
		     
			 perl check_lnL.pl > X

              filename="X"
              cat $filename | while read line
              do
                if [ "$line" = "1" ]
                 then
                 echo " "
                 echo "Higher lnL value than previous"
                 #cp highestL highestL_"$j"
			     mv -f rst ../../../04_output/tmp_result/5_alncds_rst_bmlo/rst
			     mv -f mlf ../../../04_output/tmp_result/mlf_bin"$i"
                fi;
                if [ "$line" = "0" ]
                 then
                 rm rst
                 rm mlf
                 echo " "
                 echo "Lower lnL value than previous"
                fi;
             done

			rm X
			rm $b
			rm 2base.t
			rm in.basemlg
			rm lnf
			rm rates			
			rm rst1
			rm rub			 
			}
			rm highestL
			# BASEML end #
			cd ../../../00_bin
			#-----------------------------------------------------------
			# Proc5b: BASEML output .rst -> filtered rst
			#---------------------------------------
			b="temp_b" # batchfile for each input file
			echo "../02_ctlfiles/sample_ctl/5b_34_rstfilter.ctl" >$b
			./00_mfa_convert_4_LLVMgcc_120719 < $b
			rm -r ../20.output
			rm $b
			
			#rm -rf /Volumes/data_500/runBaseml/04_output/mss_nst25/5_mstyeo_alncds_f1_rst_bmlo
		fi;
	fi;
	

	# remove temp folders #
	if [ $rmopt -eq 1 ]; then
		rm -r ../04_output/tmp_result/3_alncds_PAML
		rm -r ../04_output/tmp_result/4_alncds_PAML_3pos
		rm -r ../04_output/tmp_result/5_alncds_rst_bmlo
#		rm -r ../04_output/tmp_result/5_alncds_filter_rst 
	fi;
	
	cd ../04_output
	mv tmp_result sample_result/result_bin_"$i"
	
done

rm -r ../10.input
#----------------------------------------------------------------------------------------


