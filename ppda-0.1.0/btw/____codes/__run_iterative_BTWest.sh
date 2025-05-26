echo "start to run iterative BTWest";

# Define function called "abort".
# If this is called, then this process will be killed with a given message. 
function abort
{
   echo "$@" 1>&2
   exit 1
}

#make output directory to store the result ot BTW analysis

mkdir ../iterative_BTWest

#move input files to "____code" directory

cp ../BASEML_result/baseml_output/04_output/sample_result/result_bin_1/2c_alncds_cc_slct.MFA/concate_seq.mfa collapse_seq.txt;
cp ../BASEML_result/baseml_output/04_output/sample_result/result_bin_1/5_alncds_filter_rst/fltrst fltrst;
cp ../_concatenated_original_alignment_folder/concatenated_original_alignment.txt original_aln.txt;


# parameter setting for iterative BTW
# [HY] Load variables from control file.
control_file=$1
. $control_file

# Check if output files do not exist yet.
btw_output_file=../iterative_BTWest/anc_site_probs_all_node.txt
[ -f $btw_output_file ] && abort "$btw_output_file exist already!!!"
# TODO: Checks for other tempolary files can be added.

# Show input parameters on standard output
echo "\n[INPUT: Run iterative BTWest]"
echo "total_seq_num_in_collapse: $total_seq_num_in_collapse" # Number of sequences in the collapse sequence alignments
echo "internal_node_num: $internal_node_num" # Number of internal nodes in the input tree of BASEML
echo "poly_sp_num: $poly_sp_num" # Number of species with polymorphism samples
echo "poly_sp1_sample_start: $poly_sp1_sample_start" # Order of the first sequence of the first species in the original sequence
echo "poly_sp1_sample_end: $poly_sp1_sample_end" # Order of the last sequence of the first species in the original sequence
echo "poly_sp1_collapse1: $poly_sp1_collapse1" # Order of the first sequence of the first species in the collapse sequence
echo "poly_sp1_collapse2: $poly_sp1_collapse2" # Order of the last sequence of the first species in the collapse sequence
echo "poly_sp2_sample_start: $poly_sp2_sample_start" # Order of the first sequence of the first species in the original sequence
echo "poly_sp2_sample_end: $poly_sp2_sample_end" # Order of the last sequence of the first species in the original sequence
echo "poly_sp2_collapse1: $poly_sp2_collapse1" # Order of the first sequence of the first species in the collapse sequence
echo "poly_sp2_collapse2: $poly_sp2_collapse2" # Order of the last sequence of the first species in the collapse sequence
echo "iteration_num: $iteration_num"
echo "mutation_category_num: $mutation_category_num\n" # Number of mutation categories

# Assign polymorphism info to arrays, in which the first element represents 
# information for the first species and so on.
# This is still the same as "hard-coded" variable assignment.
# TODO: Change to recieve these arrays from control file.
declare -a poly_sample_start=($poly_sp1_sample_start $poly_sp2_sample_start)
declare -a poly_sample_end=($poly_sp1_sample_end $poly_sp2_sample_end)
declare -a poly_collapse_start=($poly_sp1_collapse1 $poly_sp2_collapse1)
declare -a poly_collapse_end=($poly_sp1_collapse2 $poly_sp2_collapse2)

mut_TC=1;
mut_TA=2;
mut_TG=3;
mut_CT=4;
mut_CA=5;
mut_CG=6;
mut_AT=7;
mut_AC=8;
mut_AG=9;
mut_GT=10;
mut_GC=11;
mut_GA=12;

#parameter setting end

#first, run BTW without weighting to make fltrst shows each site inference result (this is to make the following processes faster)
#get the position of the traget ancestral node
perl 21_read_fltrst_and_get_the_target_ancestral_node.pl $total_seq_num_in_collapse $poly_sp1_collapse1 $poly_sp1_collapse2;
perl 21_make_ancestral_site_probability_BTWnothing_method_12_mutation_categories_output_states_at_ancestral_nodes_for_new_BTW.pl $total_seq_num_in_collapse $internal_node_num $poly_sp1_sample_start $poly_sp1_sample_end $poly_sp1_collapse1 $poly_sp1_collapse2;
perl 21_recreate_fltrst_from_anc_site_probs_all_node.pl;

rm anc_site_probs_all_node.txt;
rm target_node;

#weighting starts
#replicate for the number of species with polymorphism data

for((X=0;X<=poly_sp_num-1;X++))

do

#define parameter to use

d=${poly_sample_start[$X]}
e=${poly_sample_end[$X]}
f=${poly_collapse_start[$X]}
g=${poly_collapse_end[$X]}

echo $d $e $f $g

#iterative BTW for the polymorphism of the Xth species

    #get the position of the traget ancestral node
    perl 21_read_fltrst_and_get_the_target_ancestral_node.pl $total_seq_num_in_collapse $f $g;

    #single BTWne
    perl 21_make_ancestral_site_probability_BTWne_method_12_mutation_categories_output_states_at_ancestral_nodes_for_new_BTW.pl $total_seq_num_in_collapse $internal_node_num $d $e $f $g;
    perl 21_make_estimated_SFS_BTW_method.pl $d $e;

    rm anc_site_probs.txt;
    rm anc_site_probs_all_node.txt;
    
    #iteration of BTWest
    Y=1;
    
    while [ $Y -ne $iteration_num ]

    do
         
         #run BTWest using the estimated SFS
         perl 21_make_ancestral_site_probability_BTWest_method_12_mutation_categories_output_states_at_ancestral_nodes_for_new_BTW.pl $total_seq_num_in_collapse $internal_node_num $mutation_category_num $d $e $mut_TC $mut_TA $mut_TG $mut_CT $mut_CA $mut_CG $mut_AT $mut_AC $mut_AG $mut_GT $mut_GC $mut_GA $f $g;
         rm estimated_frequency_spectrum.txt;
         #make new SFS for the next iteration
         perl 21_make_estimated_SFS_BTW_method.pl $d $e;
 
         rm anc_site_probs.txt;
        rm anc_site_probs_all_node.txt;
 
        let Y++;

    done
    
    #the last round of BTWest
    perl 21_make_ancestral_site_probability_BTWest_method_12_mutation_categories_output_states_at_ancestral_nodes_for_new_BTW.pl $total_seq_num_in_collapse $internal_node_num $mutation_category_num $d $e $mut_TC $mut_TA $mut_TG $mut_CT $mut_CA $mut_CG $mut_AT $mut_AC $mut_AG $mut_GT $mut_GC $mut_GA $f $g;
     rm estimated_frequency_spectrum.txt;
     
     rm fltrst_re;

    #create fltrst from the result of BTW for the Xth species, which will be the input of the BTW for the X+1th species
    perl 21_recreate_fltrst_from_anc_site_probs_all_node.pl;
    
    mv anc_site_probs_all_node.txt ../iterative_BTWest/anc_site_probs_all_node_"$X".txt;
    rm anc_site_probs.txt;
    
    rm target_node;

done

rm original_aln.txt;
rm collapse_seq.txt;
rm fltrst_re;
rm fltrst

#combine BTW results for the all species
mv ../iterative_BTWest/anc_site_probs_all_node_0.txt ../iterative_BTWest/anc_site_probs_all_node.txt;
perl 21_make_combine_anc_site_probs_files.pl $poly_sp_num;

rm anc_site_probs.txt;

echo "finish to run iterative BTWest";



