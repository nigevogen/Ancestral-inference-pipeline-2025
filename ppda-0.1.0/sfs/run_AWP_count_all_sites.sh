
# This script is for single allele analysis not with polymorphisms.

usage_exit() {
        echo "Usage: $0 [-p position_file] [-q position_type] [-c concat_aln_file] [-b bin_list_file] [-t rst_file] [-m mlf_file] [-o output_dir] [-i internal_node_num] [-n bootstrap_num]" 1>&2
        exit 1
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -p|--position_file)
    position_file="$2"
    shift # past argument
    shift # past value
    ;;
    -q|--position_type)
    position_type="$2"
    shift # past argument
    shift # past value
    ;;
    -c|concat_aln_file)
    concat_aln_file="$2"
    shift # past argument
    shift # past value
    ;;
    -b|---bin_list_file)
    bin_list_file="$2"
    shift # past argument
    shift # past value
    ;;
    -t|---rst_file)
    rst_file="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--mlf_file)
    mlf_file="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output_dir)
    output_dir="$2"
    shift # past argument
    shift # past argument
    ;;
    -i|--internal_node_num)
    internal_node_num="$2"
    shift # past argument
    shift # past argument
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

mkdir $output_dir
anc_prob_file=$output_dir/01_anc_prob_list.txt
joint_prob_file=$output_dir/02_anc_joint_prob_list.txt
bootstrap_results_dir=$output_dir/bootstrap_results
control_file=$output_dir/0.arguments.txt

mkdir $bootstrap_results_dir

echo "\n[ARGUMENTS]"
echo "position_file    = ${position_file}"
echo "position_type    = ${position_type}"
echo "concat_aln_file   = ${concat_aln_file}"
echo "bin_list_file     = ${bin_list_file}"
echo "rst_file          = ${rst_file}"
echo "mlf_file          = ${mlf_file}"
echo "output_dir        = ${output_dir}"
echo "internal_node_num = ${internal_node_num}"
echo "bootstrap_num     = ${bootstrap_num}"


package_dir=/Users/haruka/Documents/01_myPackages/shell/ppda

echo "from ${package_dir}/sfs/run_AWP_count_all_sites.sh\n" > $control_file
echo "[ARGUMENTS]" >> $control_file
echo "position_file    = ${position_file}" >> $control_file
echo "position_type    = ${position_type}" >> $control_file
echo "concat_aln_file   = ${concat_aln_file}" >> $control_file
echo "bin_list_file     = ${bin_list_file}" >> $control_file
echo "rst_file          = ${rst_file}" >> $control_file
echo "mlf_file          = ${mlf_file}" >> $control_file
echo "output_dir        = ${output_dir}" >> $control_file
echo "internal_node_num = ${internal_node_num}" >> $control_file
echo "bootstrap_num     = ${bootstrap_num}" >> $control_file

python3 $package_dir/btw/make_anc_prob_list_from_rst.py $rst_file $concat_aln_file $anc_prob_file

python3 $package_dir/sfs/anc_prob_list_to_joint_prob.py \
-a $anc_prob_file \
-m $mlf_file \
-j $joint_prob_file \
-i $internal_node_num

python3 $package_dir/sfs/count_mutations_by_AWP_method.py \
-o $bootstrap_results_dir \
-j $joint_prob_file \
-p $position_file \
-q $position_type \
-i $bin_list_file \
-n $bootstrap_num
