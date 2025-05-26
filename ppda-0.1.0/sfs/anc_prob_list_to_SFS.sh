
ppda_output_dir=$1

pos_table="$ppda_output_dir/1_alignment_processing/2_filtered_without_markers.csv"
bin_list_file="$ppda_output_dir/1_alignment_processing/2_filtered_without_markers/0.itemlist"
mlf_file="$ppda_output_dir/2_BTW/BASEML_result/baseml_output/04_output/sample_result/result_bin_1/mlf_bin1"
anc_prob_file="$ppda_output_dir/2_BTW/iterative_BTWest/anc_site_probs_all_node.txt"
joint_prob_file="$ppda_output_dir/3_SFS/joint_probs.txt"
summary_file="$ppda_output_dir/3_SFS/anc_state_summary.txt"
sfs_out_dir="$ppda_output_dir/3_SFS/00_all_gene_all_sites"

# TODO: These arguments should be read from control file or command line.
poly_info="RG:RG_collapse_:3,4:9:0:14..MD:MD_collapse_:5,6:10:14:21"
internal_node_num=4
bootstrap_num=1000

package_dir=/Users/haruka/Documents/01_myPackages/shell/ppda

# Compute joint probability from ancestor state probability for each position
python3 $package_dir/sfs/anc_prob_list_to_joint_prob.py \
-a $anc_prob_file \
-m $mlf_file \
-j $joint_prob_file \
-p $poly_info \
-i $internal_node_num

# Output ancestor inference summary
python3 $package_dir/sfs/output_AI_result_summary.py \
-o $summary_file \
-a $anc_prob_file \
-m $mlf_file \
-j $joint_prob_file \
-p $poly_info \
-i $internal_node_num

# Make SFS with bootstrap resampling
python3 $package_dir/sfs/count_mutations_by_AWP_method.py \
-o $sfs_out_dir \
-j $joint_prob_file \
-i $bin_list_file \
-p $pos_table \
-q table \
-n $bootstrap_num
