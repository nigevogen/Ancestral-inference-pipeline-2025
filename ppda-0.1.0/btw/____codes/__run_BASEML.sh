#path export, not necessary if baseml is in the directory already path export is done
# export PATH=$PATH:/Users/tmatsumoto/bin/;
# HY commented out here. I think it is not very good idea to modify PATH everytime running a pipeline.
# TODO: It may be better to provide baseml executable in this pipeline.

echo "start to run BASEML"; 

#make output directory
mkdir ../BASEML_result;

#copy sequences to the input folder of BASEML
cp -r -f ../_collapse_seqs_folder/ ../_AI_for_collapse/2_BASEML_and_HM/03_seq_dat/sample_data/sample_sequences/;
#copy tree file to the input folder of BASEML
cp -r -f ../_tree_folder/ ../_AI_for_collapse/2_BASEML_and_HM/03_seq_dat/;
#move bin_dat_1_index_lists to the input folder of BASEML
mv -f ../bin_dat_1_index_lists  ../_AI_for_collapse/2_BASEML_and_HM/03_seq_dat/sample_data/sample_concat_list_codons_nucleotides/;

cd ../_AI_for_collapse/2_BASEML_and_HM/01_batchfiles/batch_A/batch_B;

#run BASEML pipeline 
./sample_MFA2rst_nucleotides.sh;
cd ../../../../../;
 
rm -r BASEML_result/baseml_output;
mkdir BASEML_result/baseml_output;

#run BASEML outputs to the output directory
mv -f _AI_for_collapse/2_BASEML_and_HM/04_output/ BASEML_result/baseml_output/;

cd ____codes;

echo "finish to run BASEML"; 