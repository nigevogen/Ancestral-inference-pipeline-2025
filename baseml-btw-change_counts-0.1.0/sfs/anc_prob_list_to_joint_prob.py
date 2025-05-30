import os
import re
import argparse
import configparser
import numpy as np
import pandas as pd
from datetime import datetime
from collections import defaultdict, Counter, namedtuple

def main(
        joint_prob_file, 
        anc_prob_file, 
        mlf_file,
        internal_node_num,
        poly_info
        ):
    """ Compute joint probabilities for each mutation for every sites. Allele 
    frequency in a population is added if a mutation is segregating in a 
    population. 
    """

    start = datetime.now().isoformat()
    out_ctl_file = joint_prob_file+'.run_info.txt'
    to_control_file(
        out_ctl_file, start, 
        **{
            'joint_prob_file': joint_prob_file, 
            'anc_prob_file': anc_prob_file, 
            'mlf_file': mlf_file,
            'internal_node_num': internal_node_num, 
            'poly_info': poly_info
        }
    )
    
    calculate_joint_probability(
        joint_prob_file, anc_prob_file, mlf_file, internal_node_num, poly_info
    )

def calculate_joint_probability(
        output_file, input_path, mlf_file, internal_node_num, poly_info):
    """ Count the number of mutations in lineage-specific manner. The function
    reads an ancestral state probability file and output tab-separated list of mutation
    counts for each site. 

    Parameters
    ----------
    output_file: str
        A path to output file.
    input_path: str
        A path to ancestral state probability file output from BTW 
        ("anc_site_probs_all_node_2nd_BTW_1.txt" usually).
    mlf_file: str
        A path to mlf file. This is to get node identifiers assigned by BASEML.
    internal_node_num: int
        Number of internal nodes.
    poly_info: list
        A list of polymorphism information.
    """
    if os.path.isfile(output_file):
        raise FileExistsError(output_file)

    print('\nRead ancestor prob file', end=' ')
    poly_extant_nodes_list = [
        node for poly in poly_info for node in poly.extant_nodes]
    anc_probs_d = read_anc_site_probs(
        input_path, len(poly_info), poly_extant_nodes_list, internal_node_num)
    print('.. done')

    print('\nCount change', end=' ')
    branches = get_branch_line_from_mlf(mlf_file)
    poly_allele_num_d = {
        poly.ancestor_node: poly.allele_count for poly in poly_info
    }
    joint_probs_d = compute_joint_probs_of_all_pos(
        anc_probs_d, poly_allele_num_d, branches)
    print('.. done.')

    save_as_file(joint_probs_d, output_file)
    print(f'\nSaved to {output_file}')

def read_anc_site_probs(file_path, poly_sp_num, poly_extant_nodes_list, internal_node_num):
    """ Reads an output file of BTW. """
    d = {}
    with open(file_path, 'r') as f:
        for l in f:
            pos, poly_allele_freqs, extant, anc_probs = read_ancestor_state(
                l[:-1], poly_sp_num, poly_extant_nodes_list, internal_node_num)
            d[pos] = (poly_allele_freqs, extant, anc_probs)
            
    assert len(d) == pos+1
    return d

# TODO: Make namedtuple object for ancestral state of each pos?
def read_ancestor_state(s, poly_sp_num, poly_extant_nodes_list, internal_node_num):
    """ Reads a line for internal node nucleotide configurations and probability."""
    parts = s.split()
    pos = int(parts[0])

    # Read polymorphism allele frequency
    # poly_sp_num = len(poly_info)
    poly_info_end_ix = int(poly_sp_num*2)+1
    poly_allele_freqs = [int(f) for f in parts[1:poly_info_end_ix]]
    poly_freq_d = dict(zip(poly_extant_nodes_list, poly_allele_freqs))
    # Example: {3: 14, 4: 0, 5: 20, 6: 0}

    extant = parts[poly_info_end_ix][:-1]
    anc_list = parts[poly_info_end_ix+1:]
    n = internal_node_num+1

    # Make nested list
    anc_probs = [
        format_anc_prob(anc_list[n*i:n*i+n], internal_node_num)
            for i in range(int(len(anc_list) / n))
    ]
    
    return pos, poly_freq_d, extant, anc_probs

def format_anc_prob(anc_state, internal_node_num):
    """ Format ['G', 'A', 'G', 'G', '0.93993993993994'] to 
    ('GAGG', 0.93993993993994).
    """
    anc_state_str = ''.join(anc_state[:internal_node_num])
    prob = float(anc_state[internal_node_num])
    return (anc_state_str, prob)

def compute_joint_probs_of_all_pos(anc_probs_d, poly_allele_num_d, branches):
    """ Compute joint probabilities for each mutation type on every branches 
    across sites.
    """
    return {
        p: compute_joint_probs(
            item[0], item[1], item[2], poly_allele_num_d, branches) 
        for p, item in anc_probs_d.items()
    }

def compute_joint_probs(poly_freqs, extant, anc_probs, poly_allele_num_d, branches):
    """ Compute the sum of probabilities across internal node nucleotide 
    configurations for each mutation type on every branch on a given site. """
    scenarios = [
        (extant+anc_prob[0], anc_prob[1])
        for anc_prob in anc_probs
    ]
    tmp_d = defaultdict(list)

    # For each branch
    for anc_node, der_node in branches:
        for s in scenarios:
            anc_state, der_state = s[0][anc_node-1], s[0][der_node-1]
            prob = s[1]
            
            if anc_state == der_state:
                continue
            if prob == 0:
                continue

            # If the input anc pro file included polymorphism frequency
            if len(poly_freqs) > 0:
                # If current branch is one of polymorphism substitution lineage,
                if der_node in poly_allele_num_d:
                    # poly_allele_num_d = {9: 14, 10: 21} for mpspye tree.
                    # Key is ancestral node of polymorphism and value is total 
                    # number of alleles for the population.
                    freq = poly_allele_num_d[der_node]
                # If current branch is polymorphism lineage,
                elif der_node in poly_freqs:
                    freq = poly_freqs[der_node]
                else:
                    freq = -9
            else:
                freq = -9
            
            tmp_d[(anc_node, der_node, anc_state, der_state, freq)].append(prob)

    return {k: sum(v) for k, v in tmp_d.items()}

def get_branch_line_from_mlf(mlf_file):
    matches = []
    with open(mlf_file, 'r') as f:
        for l in f:
            if re.match(r'^\s+(\d+\.\.\d+\s+)+$', l):
               matches.append(l.rstrip()) 
    assert len(matches) == 1
    match = matches[0]

    return [[int(c) for c in a.split('..')] for a in match.split()]

def save_as_file(changes_d, out_path):
    """ Output joint probabilities to a text file. """
    lines = []
    for pos, prob_d in changes_d.items():
        if len(prob_d) == 0:
            continue
        
        # Connect allele freq with probability
        # Extract allele frequency for a given population at a position
        # freq_d[pos][ancestral_node] -> Counter({'T': 4, 'A': 10})
        # Here, ancestral_node is for a population sample
        mutations = []
        for info, prob in prob_d.items():
            anc_node, der_node, anc_state, der_state, freq = info
            mutations.append(
                f'{anc_node}..{der_node}:{anc_state}{der_state}:{freq}:{prob}'
            )
        
        lines.append('>{}\n{}'.format(pos, "\t".join(mutations)))
        
    with open(out_path, 'w') as f:
        print('\n'.join(lines), file=f)

def to_control_file(out_ctl_file, start_date, **kwargs):
    lines = [f'{k} = {v}' for k, v in kwargs.items()]

    with open(out_ctl_file, 'w') as f:
        print(f'Run date: {start_date}.', file=f)
        print('Listed parameters are used for the run.\n', file=f)
        print('\n'.join(lines), file=f)

PopulationSampleInfo = namedtuple(
    'PopulationSampleInfo', [
        'sample_prefix', # Prefix in sample seq
        'collapse_prefix', # Prefix in collapse seq
        'extant_nodes', # Order among collapse sequneces
        'ancestor_node', # Ancestral node id in a tree in mlf file (BASEML output)
        'sample_seq_start', # Position of the first allele of the population samples
        'allele_count' # Number of the population samples
    ]
)

def poly_info_parser_from_terminal(poly_info_str):
    """ Information for a population can be listed like, 
    {sample_prefix}:{prefix_in_collapse}:{order_in_collapse}:{anc_node}:
    {sample_seq_start}:{allele_count}
    If there are more than one populations, population info can be listed by 
    separating different population info with ".." mark.
    """
    if poly_info_str == '':
        return []

    populations = [l.strip() for l in poly_info_str.split('..')]
    poly_info_list = []
    for pop in populations:
        parts = pop.split(':')
        poly_info_list.append(
            PopulationSampleInfo(
                parts[0],
                parts[1],
                [int(a) for a in parts[2].split(',')],
                int(parts[3]), 
                int(parts[4]), 
                int(parts[5])
            )
        )

    return poly_info_list

def poly_info_parser_from_ctl(poly_info_str):
    if poly_info_str == '':
        return []
        
    populations = [
        l.strip() for l in poly_info_str.split('\n') if l.strip() != '']

    poly_info_list = []
    for pop in populations:
        parts = pop.split('\t')
        try:
            poly_info_list.append(
                PopulationSampleInfo(
                    parts[0], [int(a) for a in parts[1].split(',')],
                    int(parts[2]), int(parts[3]), int(parts[4])
                )
            )
        except:
            raise Exception(poly_info_str)

    return poly_info_list

def read_control_file(control_file):
    # Initialize ConfigParser object
    config = configparser.ConfigParser(
        strict=True,
        comment_prefixes=('#', ';'),
        inline_comment_prefixes=('#', ';')
    )

    # Parse control file
    paths = config.read(control_file)

    # Check number of read control files.
    if len(paths) == 0:
        raise FileNotFoundError(
        f'Specified control file, {control_file}, is not found.')
    elif len(paths) > 1:
        raise TypeError(f'Iterable {type(control_file)} is given as a control '\
            'file. Only one control file is supported.')
    assert 'SFS' in config.sections(), \
        'SFS section is not found.'

    return config['SFS']

if __name__ == '__main__':
    desc = f'Compute joint probability of each mutation on branch site '\
        'for every site.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-j", "--joint_prob_file", 
        help="A path to an output file including a list of joint "\
             "probabilities of each mutation across sites."
    )
    parser.add_argument(
        "-a", "--anc_prob_file", 
        help="A path to ancestral probability list file output from BTW."
    )
    parser.add_argument(
        "-m", "--mlf_file", 
        help="mlf file output from BASEML to get branch info."
    )
    parser.add_argument(
        "-i", "--internal_node_num", 
        help="Number of internal nodes in a sample tree used in BASEML run.",
        type=int
    )
    parser.add_argument(
        "-p", "--poly_info", 
        help="Information for population samples. "\
            "Example: RG:RG_collapse_:3,4:9:0:14..MD:MD_collapse_:5,6:10:14:21",
        nargs='?', const='', default=''
    )
    parser.add_argument(
        "-c", "--ctl_file",
        help="A path to control file. All parameters can be given via a control "\
             "file. If a control file is given, other arguments will be ignored."
    )
    args = parser.parse_args()
    if args.ctl_file:
        params = read_control_file(args.ctl_file)
        main(
            params['joint_prob_file'], 
            params['anc_prob_file'], 
            params['mlf_file'],
            int(params['internal_node_num']), 
            poly_info_parser_from_ctl(params['poly_info'])
        )

    else:
        main(
            args.joint_prob_file, 
            args.anc_prob_file,
            args.mlf_file, 
            args.internal_node_num, 
            poly_info_parser_from_terminal(args.poly_info), 
        )
