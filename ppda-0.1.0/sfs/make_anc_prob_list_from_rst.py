
import os
import numpy as np
from collections import namedtuple

AncestorState = namedtuple(
    'AncestorState',
    ['index', 'obs_cnt', 'extant_config', 'anc_probs', 'total_prob']
)

def main(rst_file, aln_file, out_file):
    if os.path.isfile(out_file):
        raise FileExistsError(out_file)

    anc_state_d = read_rst_file(rst_file)
    seqnames, seq_matrix = read_concatenate_aln_file(aln_file)

    anc_state_list = map_anc_prob_to_every_site(anc_state_d, seq_matrix)
    out_lines = [
        anc_state_to_str(anc, i, 'X', 'X') 
        for i, anc in enumerate(anc_state_list)]

    with open(out_file, 'w') as f:
        print('\n'.join(out_lines), file=f)

def read_rst_file(rst_file):
    anc_state_d = {}
    start_parse = False
    check_site_num = False
    obs_site_cnt = 0

    with open(rst_file, 'r') as f:
        for l in f:
            if l.strip() == '':
                continue

            if l.startswith('Pattern Freq   Data:'):
                start_parse = True

            elif l.startswith('List of extant and reconstructed sequences'):
                start_parse = False
                check_site_num = True

            elif start_parse:
                anc_state = parse_anc_state_str(l.rstrip())
                obs_site_cnt += anc_state.obs_cnt
                anc_state_d[anc_state.extant_config] = anc_state

            elif check_site_num:
                exp_site_cnt = int(l.strip().split()[1])
                if obs_site_cnt != exp_site_cnt:
                    raise ValueError(
                    f'Wrong total site count: {obs_site_cnt} != {exp_site_cnt}')

                check_site_num = False

    return anc_state_d

def parse_anc_state_str(anc_state_str):
    parts = anc_state_str.split(':')
    keys = parts[0].split()
    try:
        anc_states = parts[1].split('(total')[0].split()
    except IndexError:
        raise IndexError(anc_state_str)

    anc_state_d = dict(zip(
        [a for i, a in enumerate(anc_states) if i % 2 == 0], 
        [float(a[1:-1]) for i, a in enumerate(anc_states) if i % 2 == 1]
    ))
    total_probs = anc_state_str.split(':')[1].split('(total')[1].split()
    assert len(total_probs) == 1
    total_prob = float(total_probs[0][:-1])
    obs_total_prob = sum(anc_state_d.values())

    if round_num(obs_total_prob, 1) != round_num(total_prob, 1):
        msg = f'Wrong total probability: '\
            f'{obs_total_prob} != {total_prob} at '\
            f'line "{anc_state_str}".'
        raise ValueError(msg)

    return AncestorState(
        int(keys[0]), 
        int(keys[1]),
        keys[2],
        anc_state_d,
        total_prob
    )

def round_num(a, ndigits=2):
    n = 10 ** ndigits
    return (a * n * 2 + 1) // 2 / n

def read_concatenate_aln_file(aln_file):
    seqnames = []
    seqs = []

    with open(aln_file, 'r') as f:
        for l in f:
            if l.startswith('>'):
                seqnames.append(l[1:-1])
            else:
                seqs.append(list(l[:-1]))

    return seqnames, np.array(seqs)

def map_anc_prob_to_every_site(anc_state_d, seq_matrix):
    out_list = [
        anc_state_d[''.join(config)]
        for config in seq_matrix.T
    ]
    return out_list

def anc_state_to_str(anc_state, site_index, allele_freq1, allele_freq2):
    items = [
        str(site_index),
        allele_freq1, 
        allele_freq2,
        anc_state.extant_config+':', 
        '\t'.join([
            '\t'.join(list(anc_states)+[str(prob)]) 
            for anc_states, prob in anc_state.anc_probs.items()
        ])
    ]
    return '\t'.join(items)

if __name__ == '__main__':
    import sys
    rst_file = sys.argv[1]
    aln_file = sys.argv[2]
    out_file = sys.argv[3]
    print('\nCAUTION: Please make sure that the given alignment file is '\
          'the one input to BASEML')
    print('\nINPUT')
    print(f'rst file: {rst_file}')
    print(f'aln file: {aln_file}')
    print(f'output file: {out_file}')

    main(rst_file, aln_file, out_file)

