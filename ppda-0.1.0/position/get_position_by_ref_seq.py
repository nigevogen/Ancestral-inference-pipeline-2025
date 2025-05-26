import os
import numpy as np
import argparse
from collections import namedtuple

def main(output_pos_file, unfilt_aln_dir, ref_seq_prefix, 
         first_cod_pos, last_cod_pos, min_cod_num=None, max_cod_num=None):

    out_pos_d = collect_output_positions(
        unfilt_aln_dir, ref_seq_prefix, first_cod_pos, last_cod_pos,
        min_cod_num, max_cod_num)
    
    lines = [
        '>{}\n{}'.format(k, '\n'.join([str(a) for a in v]))
        for k, v in out_pos_d.items()
    ]
    with open(output_pos_file, 'w') as f:
        print('itemnum:', len(lines), file=f)
        print('\n'.join(lines), file=f)

def collect_output_positions(
        unfilt_aln_dir, ref_seq_prefix, first_cod_pos, last_cod_pos, 
        min_cod_num=None, max_cod_num=None):
    items = read_item_list(os.path.join(unfilt_aln_dir, '0.itemlist'), int)
    out_pos_d = {}

    for item in items.item_list:
        fpath = os.path.join(unfilt_aln_dir, items.prefix+str(item)+items.suffix)
        out_pos = get_ref_seq_positions_from_aln_file(
            fpath, ref_seq_prefix, first_cod_pos, last_cod_pos, 
            min_cod_num, max_cod_num)
        if len(out_pos) > 0:
            out_pos_d[item] = out_pos

    return out_pos_d

def read_aln_file_to_matrix(aln_file):
    seqnames = []
    seqs = []

    with open(aln_file, 'r') as f:
        for l in f:
            if l.startswith('>'):
                seqnames.append(l[1:-1])
            else:
                seqs.append(list(l[:-1]))

    return seqnames, np.array(seqs)

def get_ref_seq_positions(
        ref_seq, first_cod_pos, last_cod_pos, 
        min_cod_num=None, max_cod_num=None):
    rm_gap = [a for a in ref_seq if a != '-']
    if min_cod_num != None:
        if len(rm_gap) // 3 < min_cod_num:
            return []

    if max_cod_num != None:
        if len(rm_gap) // 3 > max_cod_num:
            return []

    output_pos = []
    ref_pos = 0

    for aln_pos in range(0, len(ref_seq), 3):
        ref_cod = ''.join(ref_seq[aln_pos:aln_pos+3])

        if '-' not in ref_cod:        
            if (first_cod_pos-1)*3 <= ref_pos and ref_pos < last_cod_pos*3:
                output_pos.append(aln_pos+2)
                
            ref_pos += 3

    return output_pos

def get_ref_seq(seqnames, ref_seq_prefix, seq_matrix):
    ref_seq_ixs = [i for i, a in enumerate(seqnames) if a.startswith(ref_seq_prefix)]
    assert len(ref_seq_ixs) == 1
    ref_seq_ix = ref_seq_ixs[0]
    
    return seq_matrix[ref_seq_ix]

def get_ref_seq_positions_from_aln_file(
        aln_path, ref_seq_prefix, first_cod_pos, last_cod_pos,
        min_cod_num=None, max_cod_num=None):
    
    seq_names, seq_matrix = read_aln_file_to_matrix(aln_path)
    ref_seq = get_ref_seq(seq_names, ref_seq_prefix, seq_matrix)
    
    return get_ref_seq_positions(
        ref_seq, first_cod_pos, last_cod_pos, min_cod_num, max_cod_num)

ItemList = namedtuple(
    'ItemList', 
    ['name', 'prefix', 'suffix', 'item_list']
)
def read_item_list(
        item_list_path:str,
        apply_func=None
        ):
    if apply_func == None:
        apply_func = lambda s: s

    item_list = []
    with open(item_list_path, 'r') as f:
        for i, l in enumerate(f):
            l = l.strip()
            if i == 0:
                list_num = int(l)
                assert list_num == 1, \
                    f'1 is expected but {list_num} is given as list num.'
                continue
            if l.startswith('prefix:'):
                prefix = l.split('prefix:')[1].strip()
                continue

            if l.startswith('suffix'):
                suffix = l.split('suffix:')[1].strip()
                continue

            if l.startswith('>'):
                list_name = l[1:].split('\t')[0]
                itemnum = int(l.split('\t')[1])
                continue

            else:
                item_list.append(apply_func(l))
    
    assert len(item_list) == itemnum
    return ItemList(list_name, prefix, suffix, item_list)

if __name__ == '__main__':
    desc = 'Output basic statistics for all sites used for ancestor inference.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-o", "--output_pos_file", 
        help="A path to output position list file. This can be used for "\
             "filtering position map file."
    )
    parser.add_argument(
        "-u", "--unfilt_aln_dir", 
        help="A path to unfiltered alignment directory."
    )
    parser.add_argument(
        "-r", "--ref_seq_prefix", 
        help="Prefix of reference sequence."
    )
    parser.add_argument(
        "-f", "--first_cod_pos", 
        help="Start position of codon to include in output. If 1 is given, "\
             "codons from the first codon ('ATG' in many cases) will be output.",
        type=int
    )
    parser.add_argument(
        "-l", "--last_cod_pos", 
        help="Last position of codon to include in output. If 50 is given, "\
             "codons until and including the 50th codon will be output.",
        type=int
    )
    parser.add_argument(
        "-m", "--min_cod_num", 
        help="Minimum length of CDS to output. If CDS length is less than a "\
             "given value, CDS will not be in an output list.",
        type=int
    )
    parser.add_argument(
        "-M", "--max_cod_num", 
        help="Maximum length of CDS to output. If CDS length is greater than a "\
             "given value, CDS will not be in an output list.",
        type=int
    )
    args = parser.parse_args()
    main(
        args.output_pos_file, 
        args.unfilt_aln_dir, 
        args.ref_seq_prefix, 
        args.first_cod_pos, 
        args.last_cod_pos,
        args.min_cod_num, 
        args.max_cod_num, 
    )
