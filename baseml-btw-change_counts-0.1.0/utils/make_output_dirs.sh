# Make output directories for the pipeline

# Get a file path to output directory
output_dir=$1

# TODO: Need to generalize
# Please change this value to the location of this folder on your computer. 
package_dir=/Users/<username>/<thisfolder>

# Make output directory for alignments
mkdir "$output_dir/1_alignment_processing"

# Make output directory for BTW
mkdir "$output_dir/2_BTW"
cp "$package_dir/btw/____codes.zip" "$output_dir/2_BTW"
unzip -q "$output_dir/2_BTW/____codes.zip" -d "$output_dir/2_BTW"

cp "$package_dir/btw/_AI_for_collapse.zip" "$output_dir/2_BTW"
unzip -q "$output_dir/2_BTW/_AI_for_collapse.zip" -d "$output_dir/2_BTW"

cp "$package_dir/btw/run_BTW_pipeline.sh" "$output_dir/2_BTW"

mkdir "$output_dir/2_BTW/_seq_list_to_be_concatenated"
mkdir "$output_dir/2_BTW/_seqs_folder"
mkdir "$output_dir/2_BTW/_tree_folder"

# rm "$output_dir/2_BTW/____codes.zip"
# rm "$output_dir/2_BTW/_AI_for_collapse.zip"

# Make output directory for SFS
mkdir "$output_dir/3_SFS"
