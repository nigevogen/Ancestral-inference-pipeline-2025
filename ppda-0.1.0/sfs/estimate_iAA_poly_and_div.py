""" Hard-coded script to estimate polymorphism and divergence for individual 
amino acids. """

import os
import argparse

AMINO_ACIDS = [
    'F', 'Y', 'H', 'Q', 'N', 'K', 'D', 
    'E', 'C', 's2', 'V', 'S', 'P', 'T', 'A', 'G'
]

def main(analysis_dir, pos_dir_name, ref_posmap_list_path, rep_num):
    # Make position list for iAA
    print('\n### Make position lists for iAA ###')
    script = '/Users/haruka/Documents/01_myPackages/shell/ppda/position/'\
             'get_position_by_amino_acid_class.py'
    aln_dir = os.path.join(
        analysis_dir, '1_sample_sequences/1_unfiltered_with_markers')
    pos_list_dir = os.path.join(
        analysis_dir, f'4_AWP_mutation_count/{pos_dir_name}')
    if not os.path.isdir(pos_list_dir):
        os.mkdir(pos_list_dir)
    out_prefix = '01_iAA_'
    out_suffix = '.3rd_pos.map.list'

    command = f'python3 {script} -a {aln_dir} -ap Dmel_ -as .mpspye.aln '\
        f'-m AA_marker -o {pos_list_dir} -op {out_prefix} -os {out_suffix}'
    os.system(command)
    
    os.mkdir(os.path.join(analysis_dir, f'4_AWP_mutation_count/01_iAA/'))

    for aa in AMINO_ACIDS:
        # Filter sites
        print('\n### Map to filtered alignments ###')
        pos_list_path = os.path.join(pos_list_dir, f'{out_prefix}{aa}{out_suffix}')
        out_pos_list_path = os.path.join(pos_list_dir, f'{out_prefix}{aa}{out_suffix}')
        script = '/Users/haruka/Documents/01_myPackages/shell/ppda/position/'\
                'filter_position_map_list.py'

        command = f'python3 {script} -o {out_pos_list_path} -i {ref_posmap_list_path} '\
            f'-f {pos_list_path} -m retain -n 0'
        os.system(command)

        # Bootstrap and group mutation counts
        out_dir = os.path.join(
            analysis_dir, f'4_AWP_mutation_count/01_iAA/{out_prefix}{aa}')
        os.mkdir(out_dir)
        
        joint_path = os.path.join(
            analysis_dir, '4_AWP_mutation_count/join_prob_with_freq.txt')
        bin_list_path = os.path.join(
            analysis_dir, '1_sample_sequences/5_non_collapse_concat_seq_3rd_pos/baseml_bin_list')
        unfilt_cds_dir = os.path.join(
            analysis_dir, '1_sample_sequences/1_unfiltered_with_markers')
        
        script = '/Users/haruka/Documents/01_myPackages/shell/ppda/sfs/'\
                'count_mutations_by_AWP_method.py'
        command = f'python3 {script} -o {out_dir} -j {joint_path} '\
            f'-p {out_pos_list_path} -q posmap -i {bin_list_path} '\
            f'-s {unfilt_cds_dir} -g a -n {rep_num}'

        os.system(command)

if __name__ == '__main__':
    desc = 'Estimate polymorphism and divergence for individual amino acids.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        "-a", "--analysis_dir", 
        help="A path to analysis directory."
    )
    parser.add_argument(
        "-p", "--pos_dir_name", 
        help="A name of position list directory. This has to be under "\
             "4_AWP_mutation_count."
    )
    parser.add_argument(
        "-r", "--ref_posmap_list_path", 
        help="A path to position map list for all sites."
    )
    parser.add_argument(
        "-n", "--rep_num", 
        help="A number of bootstrap sampling.", type=int
    )
    
    args = parser.parse_args()
    main(
        args.analysis_dir, 
        args.pos_dir_name, 
        args.ref_posmap_list_path, 
        args.rep_num
    )

