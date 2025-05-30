echo "start to make concatenated original alignment"; 

#make output directory
mkdir ../_concatenated_original_alignment_folder;

#run perl script
perl make_concatenated_original_seq.pl
  
#move output to the directory
mv concatenated_original_alignment.txt ../_concatenated_original_alignment_folder/;

echo "finish to make concatenated original alignment"; 