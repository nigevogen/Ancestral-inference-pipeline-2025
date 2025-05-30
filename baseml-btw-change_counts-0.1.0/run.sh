
# Define function called "abort".
# If this is called, then this process will be killed with a given message. 
function abort
{
   echo "$@" 1>&2
   exit 1
}

# Get a path to a control file from command line
control_file=$1

# Asssign a path to pipeline package
package_dir=/Users/<username>/<thisfolder> # Please change this value to the location of this folder on your computer. 

# Make tempolary control file that can be loaded from shell script.
python3 $package_dir/utils/make_template_control_file.py $control_file 0.shell_vars.ini

# Load variables declared in 0.shell_vars.ini
# TODO: this file name can be with a random number to avoid conflict from other 
# processes.
. 0.shell_vars.ini

echo "\nINPUT"
echo "Analysis directory: $ppda_output_dir"
echo "Input tree file: $sample_tree_file\n"

# Move shell control file pipeline output directory
mv 0.shell_vars.ini "$ppda_output_dir/0.shell_vars.ini"

# The following arguments are loaded via 0.shell_vars.ini,
# - ppda_output_dir
# - site_type
# - sample_tree_file
# - iteration_num
# - internal_node_num
# - mutation_category_num
# - poly_sp_num
# - poly_sp1_sample_start
# - poly_sp1_sample_end
# - poly_sp1_collapse1
# - poly_sp1_collapse2
# - poly_sp2_sample_start
# - poly_sp2_sample_end
# - poly_sp2_collapse1
# - poly_sp2_collapse2
# - single_allele_sp_num
# - single_allele_sp1_sample_pos
# - single_allele_sp2_sample_pos
# - total_seq_num_in_collapse

# Assign variables for directory names
aln_dir="$ppda_output_dir/1_alignment_processing"
btw_dir="$ppda_output_dir/2_BTW"
sfs_dir="$ppda_output_dir/3_SFS"

echo $aln_dir
echo $btw_dir
echo $sfs_dir

# Assign a variable for sample alignment directory.
# [HY] Other options may be added in the future (e.g. not only 3rd position 
# but also codon reconstruction)
case "$site_type" in
  "nucleotide" ) sample_aln_dir="$aln_dir/2_filtered_without_markers";;
  "codon" ) sample_aln_dir="$aln_dir/3_filtered_without_markers_3rd_pos";;
  * ) abort "Unknown site_type is found: $site_type. Please specify codon or nucleotide. Process is aborted." ;;
esac

# Make output directories
$package_dir/utils/make_output_dirs.sh "$ppda_output_dir"

# Make temporaly control file only for alignment filtering process
python3 $package_dir/utils/make_filtering_ctl_file_from_pipeline_ctl.py \
"$control_file" \
"$aln_dir/0.alignment.ini"

# Filter alignments
python3 $package_dir/alignment/filter_alignments.py "$aln_dir/0.alignment.ini"

# Copy sample tree file.
cp $sample_tree_file "$btw_dir/_tree_folder/sample.trees"

# Copy filtered alignment files to 2_BTW/seqs_folder and output other list files
python3 $package_dir/alignment/copy_sample_seqs_to_BTW_dir.py "$sample_aln_dir" "$btw_dir"

# Run BTW
cd $btw_dir
time ./run_BTW_pipeline.sh "$ppda_output_dir/0.shell_vars.ini"

# Make SFS
time /Users/haruka/Documents/01_myPackages/shell/ppda/sfs/anc_prob_list_to_SFS.sh "$ppda_output_dir"
