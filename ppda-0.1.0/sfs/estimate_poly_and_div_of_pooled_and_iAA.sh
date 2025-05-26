
usage_exit() {
        echo "Usage: $0 [-a analysis_dir] [-n bootstrap_num]" 1>&2
        exit 1
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -a|--analysis_dir)
    analysis_dir="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--bootstrap_num)
    bootstrap_num="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help) usage_exit
    ;;
    *) usage_exit   # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

mkdir $analysis_dir/4_AWP_mutation_count/00_all_genes_all_sites
control_file=$analysis_dir/4_AWP_mutation_count/00_all_genes_all_sites/0.arguments.txt

package_dir=/Users/haruka/Documents/01_myPackages/shell/ppda/

echo "from ${package_dir}/sfs/estimate_poly_and_div_of_pooled_and_iAA.sh\n" > $control_file
echo "[ARGUMENTS]" >> $control_file
echo "analysis_dir      = ${analysis_dir}" >> $control_file
echo "bootstrap_num     = ${bootstrap_num}" >> $control_file

# Make joint prob with frequency from BASEML output
sample_seq_file=$analysis_dir/1_sample_sequences/5_non_collapse_concat_seq_3rd_pos/concat_non-collapse.aln
mlf_file=$analysis_dir/2_baseml_result/04_output/sample_result/result_bin_1/mlf_bin1
anc_prob_file=$analysis_dir/3_btw_result/anc_site_probs_all_node_2nd_BTW_1.txt
joint_prob_file=$analysis_dir/4_AWP_mutation_count/join_prob_with_freq.txt
internal_node_num=4

python3 $package_dir/sfs/anc_prob_list_to_joint_prob.py \
-a $anc_prob_file \
-m $mlf_file \
-j $joint_prob_file \
-i $internal_node_num \
-s $sample_seq_file \
-p RG:RG_collapse_:3,4:9:0:14..MD:MD_collapse_:5,6:10:14:21

# Estimate SFS and divergence for all codons
all_codon_dir=$analysis_dir/4_AWP_mutation_count/00_all_genes_all_sites/
position_file=$analysis_dir/1_sample_sequences/1_unfiltered_with_markers.csv
bin_list_file=$analysis_dir/1_sample_sequences/5_non_collapse_concat_seq_3rd_pos/baseml_bin_list

python3 $package_dir/sfs/count_mutations_by_AWP_method.py \
-o $all_codon_dir \
-j $joint_prob_file \
-p $position_file \
-q table \
-i $bin_list_file \
-n $bootstrap_num

# Estimate SFS and divergence for iAA
ref_posmap_list_path=$all_codon_dir/1_unfiltered_with_markers.csv.posmap.list

python3 $package_dir/sfs/estimate_iAA_poly_and_div.py \
-a $analysis_dir \
-p _position_filter_list \
-r $ref_posmap_list_path \
-n $bootstrap_num
