# setting --------------------------------------------------------------------------------
run_mode=0							# 0: sample of codon analysis
                                    # 1: sample of nucleotide analysis
                                    # 2: sample of codon bootstrap analysis
                                    # 3: your analysis
#-----------------------------------------------------------------------------------------


# define output folder of the run
cd ../..

count=1
while [ -e "OUTPUT_$count" ]
do
 $((count++));
done

output_folder="OUTPUT_$count"

# if the previous run is stopped on the way, move the output to "OUTPUT_$count"
# "04_output" is always for the current run

if [ -e "04_output" ]; then
 echo "previous run was not finished"

 mv 04_output $output_folder
 
 echo "output of the previous run is moved to $output_folder"
 
 $((count++));
 output_folder="OUTPUT_$count"
fi

# start BASEML analysis
cd 01_batchfiles/batch_A/batch_B

case "$run_mode" in
		0) ./sample_MFA2rst_codons.sh
		   ./sample_rst2HM_codons.sh;;
		1) ./sample_MFA2rst_nucleotides.sh
		   ./sample_rst2HM_nucleotides.sh;;
		2) ./sample_MFA2rst_codons_bootstrap_analysis.sh
		   ./sample_rst2HM_codons_bootstrap.sh;;
		3) ./your_run_MFA2rst.sh
		   ./your_run_rst2HM.sh;;
esac

   
cd ../../..

# if the analysis finished, output is moved to "OUTPUT_$count"
mv 04_output $output_folder