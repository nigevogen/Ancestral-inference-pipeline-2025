#-----------------------------------------------------------------------------------------
#- copy the following code to pause the process
#- CODE: 
#-     echo ""; echo "Press enter to continue"; read Wait;
#-----------------------------------------------------------------------------------------
# setting --------------------------------------------------------------------------------
#PAMLver=1								# 0: HM; 1: 12*12 substitutions from EMC; 2: both;
rmopt=1									# temp folder removing option
CDSopt=1								# 0: nucleotides; 1: codons; now 0 is not available
rstopt=1								# 0: bmlo rst; 1: filtered rst; now 0 is not available
filteropt=1						# 1: allcod; 2: 2f10cDcRY; 3, 4f6cN; 4: 2f10cDcAA; 5: 4f6cAA;
subopt=(1 1 1)							# subset choice option (nh0, nh1 or nh4) 
#-----------------------------------------------------------------------------------------

cd ../../../00_bin
mkdir ../10.input

for((i=1;i<=100;i++))
do
	# pick up a bin output folder (folder name will be changed back later) #
	cd ../04_output
	mv sample_result/result_bin_"$i" tmp_result
	mkdir tmp_result/8_alncds_HM
	cd ../00_bin
	#------------------------------------------------------------------------------
	# subsets model1
	#---------------------------------------------------------------	
	if [ ${subopt[0]} -eq 1 ]; then
		#-----------------------------------------------------------
		# Proc9: BASEML output .rst or filtered rst -> HitMatrix (nh0)
		#---------------------------------------
		b="temp_b" # batchfile for each input file
		if [ $CDSopt -eq 1 ]; then
			case "$rstopt" in
				0) echo "../02_ctlfiles/sample_ctl/9_35_rst2HM.ctl" >$b;;
				1) echo "../02_ctlfiles/sample_ctl/9_35_frst2HM.ctl" >$b;;
			esac
		fi;
		if [ $CDSopt -eq 0 ]; then
			case "$rstopt" in
				0) echo "../02_ctlfiles/sample_ctl/9_38_rst2nHM.ctl" >$b;;
				1) echo "../02_ctlfiles/sample_ctl/9_38_frst2nHM.ctl" >$b;;
			esac
		fi;
		./00_mfa_convert_4_LLVMgcc_120719 < $b
		rm -r ../20.output
		rm $b
	fi;
	

	# remove temp folders #
	if [ $rmopt -eq 1 ]; then
		case "$rstopt" in
			0) rm -r ../04_output/tmp_result/5_alncds_rst_bmlo;;
			1) rm -r ../04_output/tmp_result/5_alncds_filter_rst;;
		esac
	fi;
	
	cd ../04_output
	mv tmp_result sample_result/result_bin_"$i"
	
done
rm -r ../10.input



#130830matsumoto
#----------------------------------------------------------------------------------------
#if [ $PAMLver -eq 1 -o $PAMLver -eq 2 ]; then

        #-----------------------------------------------------------
		# use output of new_PAML and calculate numbers of 12 ATGC changes
		#---------------------------------------

#cd ../../../04_output;

#for((i=1;i<6+1;i++))
#do

#    mkdir -p newPAML/nh4_pos1/bin_"$i";
#    mkdir -p newPAML/nh4_pos2/bin_"$i";
#    mkdir -p newPAML/nh4_pos3/bin_"$i";
    
	#------------------------------------------------------------------------------
	# subsets nh4  newPAML works only in nh4
	#---------------------------------------------------------------	
#	if [ ${subopt[0]} -eq 1 ]; then
		#-----------------------------------------------------------
		# Proc10: mfl -> 12ATGC changes (nh4)
		#---------------------------------------
		
#		mv mlf_nh4_pos1_$i mlf;
		
#		cd ../00_bin;
#		./130905newPAML-12changes;

		
#		cd ../04_output;
		
#		mv node* newPAML/nh4_pos1/bin_"$i";
		
#		mv mlf_nh4_pos2_$i mlf;
		
#		cd ../00_bin;
#		./130905newPAML-12changes;

		
#		cd ../04_output;
		
#		mv node* newPAML/nh4_pos2/bin_"$i";
		
#		mv mlf_nh4_pos3_$i mlf;
		
#		cd ../00_bin;
#		./130905newPAML-12changes;

		
#		cd ../04_output;
		
#		mv node* newPAML/nh4_pos3/bin_"$i";
#	fi;
	
#	cd ../04_output;
#	rm mlf;

	
#done

#fi;




