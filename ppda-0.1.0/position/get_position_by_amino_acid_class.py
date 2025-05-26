
import os
import argparse
import alignmentrs as al

AMINO_ACIDS = dict([
    (aa, aa.upper()+'2') if aa.upper() != aa else (aa, aa)
    for aa in 'FYHQNKDECsVSPTAG'
])

def main(
        aln_dir:str, in_prefix:str, in_suffix:str, 
        marker_name:str, output_dir:str, out_prefix:str, out_suffix:str):
    aln_d = read_alignments_with_AA_marker(aln_dir, in_prefix, in_suffix)
    aa_pos_d = get_positions_of_individual_amino_acid(aln_d, marker_name)
    output_as_position_list(aa_pos_d, output_dir, out_prefix, out_suffix)

def read_alignments_with_AA_marker(aln_dir, prefix, suffix):
    file_list = [l for l in os.listdir(aln_dir) if l.endswith(suffix)]
    return {
        fname.split(prefix)[1].split(suffix)[0]: \
            al.fasta_file_to_alignment(os.path.join(aln_dir, fname), fname)
        for fname in file_list
    }

def get_positions_of(seq, match_char, thirdpos=True):
    if thirdpos:
        return [
            i for i, site in enumerate(seq) 
                if site == match_char and i % 3 == 2
        ]
    
    return [i for i, site in enumerate(seq) if site == match_char]

def get_positions_of_individual_amino_acid(aln_d, marker_name):
    aa_pos_d = dict.fromkeys(AMINO_ACIDS.keys())

    for aa in AMINO_ACIDS.keys():
        tmp_d = {}

        for k, aln in aln_d.items():
            aa_seq = aln.get_samples([marker_name], match_prefix=True)\
                        .sample_sequences[0]
            tmp_d[k] = get_positions_of(aa_seq, aa)

        aa_pos_d[aa] = tmp_d

    return aa_pos_d

def output_as_position_list(aa_pos_d, output_dir, prefix, suffix):
    for aa, pos_d in aa_pos_d.items():
        out_path = os.path.join(output_dir, prefix + AMINO_ACIDS[aa] + suffix)
        
        out_lines = [
            '>{}\n{}'.format(
                name, '\n'.join([str(p) for p in pos_list])
            )
            for name, pos_list in pos_d.items()
        ]
        
        with open(out_path, 'w') as f:
            print('itemnum:', len(out_lines), file=f)
            print('\n'.join(out_lines), file=f)
            
if __name__ == '__main__':
    desc = 'Make position lists of third positions of codons for individual '\
           'amino acids.'
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        "-a", "--aln_dir", 
        help="A path to input directory for alignments with amino acid (AA) "\
             "marker seq."
    )
    parser.add_argument(
        "-ap", "--aln_prefix", 
        help="A str of prefix of input file name."
    )
    parser.add_argument(
        "-as", "--aln_suffix", 
        help="A str of suffix of input file name."
    )
    parser.add_argument(
        "-m", "--marker", 
        help="A str of marker name for AA marker sequence."
    )
    parser.add_argument(
        "-o", "--out_dir", 
        help="A path to output directory for position list files."
    )
    parser.add_argument(
        "-op", "--out_prefix", 
        help="A str of prefix of output file name."
    )
    parser.add_argument(
        "-os", "--out_suffix", 
        help="A str of suffix of output file name."
    )
    
    args = parser.parse_args()
    main(
        args.aln_dir, 
        args.aln_prefix, 
        args.aln_suffix, 
        args.marker, 
        args.out_dir, 
        args.out_prefix, 
        args.out_suffix
    )
