#parameter setting
 # [HY] Load variables
control_file=$1
. $control_file

echo "\n[INPUT: Make collapse seqs]"
echo "poly_sp_num: $poly_sp_num"
echo "poly_sp1_sample_start: $poly_sp1_sample_start"
echo "poly_sp1_sample_end: $poly_sp1_sample_end"
echo "poly_sp2_sample_start: $poly_sp2_sample_start"
echo "poly_sp2_sample_end: $poly_sp2_sample_end"
echo "single_allele_sp_num: $single_allele_sp_num"
echo "single_allele_sp1_sample_pos: $single_allele_sp1_sample_pos"
echo "single_allele_sp2_sample_pos: $single_allele_sp2_sample_pos\n"

#parameter setting end

Count=0
expectedTotal=0  
dirpath="../_seqs_folder"     # input folder (same directory as .sh file)
collapsedDirpath="../_collapse_seqs_folder"  # output folder (to be created)

mkdir "$collapsedDirpath";    # create requested output folder

echo "start to make collapse sequences"; 
while IFS='' read -r line || [[ -n "$line" ]]; do                    
    if [[ $line == itemnum\:* ]] && [ $expectedTotal == 0 ]; then    
        expectedTotal=$(echo "$line" | cut -d' ' -f 2 | tr -d '\n')
        #echo $expectedTotal;
    elif [ $expectedTotal > 0 ]; then
        if [ -e "$dirpath/$line" ]; then

                # echo "$line"; 

                cp "$dirpath/$line" data.ffn; # search file with name "$dirpath/$line" and rename it to "data.ffn", which is the input file name of the perl script below
             
                perl make_collapse_seqs_test_nucleotide.pl $poly_sp_num $poly_sp1_sample_start $poly_sp1_sample_end $poly_sp2_sample_start $poly_sp2_sample_end $single_allele_sp_num $single_allele_sp1_sample_pos $single_allele_sp2_sample_pos; # run perl script to make collapse sequences
                mv ndata.ffn "$collapsedDirpath/$line";  # rename the output file name of the perl script, ndata.ffn, to the original name shown in 0.filelist and move it to "$collapsedDirpath/"
             
                rm data.ffn;
                (( Count++ ));
        fi
    fi
# done < "$1"
done < "$dirpath/0.filelist" # filename of the sequence list is given here

echo "finish to make collapse sequences"; 
