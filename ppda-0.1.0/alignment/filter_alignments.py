import os
import re
import sys
import numpy as np
import pandas as pd
import configparser
import alignmentrs as al
from random import shuffle
from copy import deepcopy
from collections import OrderedDict, namedtuple
from warnings import warn
from datetime import datetime

STOPCODONS = {'TAG', 'TGA', 'TAA'}

BASES = 'TCAG'
CODONS = [a+b+c for b in BASES for a in BASES for c in BASES]
STOP_CODONS = ['TGA', 'TAG', 'TAA']
AMINO_ACIDS = 'ARNDCQEGHILKMFPSTWYV'

transl_code_a = 'F2L6I3M1V4S4P4T4A4Y2*2H2Q2N2K2D2E2C2*1W1R4S2R2G4'
transl_a = ''.join([str(a*int(b)) for a, b in zip(transl_code_a[::2], transl_code_a[1::2])])
GENETIC_CODE_a = OrderedDict(zip(CODONS, transl_a))

# S4 s2 are considered separate categories (Other six fold is in one category)
transl_code_b = 'F2L6I3M1V4S4P4T4A4Y2*2H2Q2N2K2D2E2C2*1W1R4s2R2G4'
transl_b = ''.join([str(a*int(b)) for a, b in zip(transl_code_b[::2], transl_code_b[1::2])])
GENETIC_CODE_b = OrderedDict(zip(CODONS, transl_b))

# ------- main function ------- #

def main(
        input_aln_dir:str, 
        input_aln_list_file:str, 
        aln_out_dir:str, 
        site_type:str,
        sample_order:list,
        marker_kw:str, 
        markers_to_use:list,
        poly_info:list,
        mt2_seg_allele_filter_mode:str, 
        only_conserved_aa:bool, 
        genetic_code_type:str,
        assert_file_exist:bool, 
        trim_kw:str,
        ):
    
    start = datetime.now().isoformat()

    assert genetic_code_type in {'a', 'b'}, \
        f'Unknown genetic code type "{genetic_code_type}" is found. '\
         'Please specify "a" or "b".'
    # Prepare output directory
    dir_paths = check_output_paths(aln_out_dir, site_type)

    # Output read control file
    out_ctl_path = os.path.join(aln_out_dir, '0.read_ctl.txt')
    to_control_file(
        out_ctl_path, 
        start,
        input_aln_dir=input_aln_dir, 
        input_aln_list_file=input_aln_list_file, 
        aln_out_dir=aln_out_dir, 
        site_type=site_type,
        sample_order=sample_order,
        marker_kw=marker_kw, 
        markers_to_use=markers_to_use,
        poly_info=poly_info, 
        mt2_seg_allele_filter_mode=mt2_seg_allele_filter_mode,
        only_conserved_aa=only_conserved_aa, 
        genetic_code_type=genetic_code_type,
        assert_file_exist=assert_file_exist, 
        trim_kw=trim_kw
    )

    # Read input alignments
    aln_d, input_items = read_alignments_from_directory(
        input_aln_dir, input_aln_list_file, marker_kw, markers_to_use,
        assert_file_exist)
    check_sample_size(aln_d)

    # Extract population sample names
    poly_name_in_collapse = [
        [poly.sample_prefix, poly.collapse_prefix] 
        for poly in poly_info
    ]
    poly_sample_prefix_list = [l[0] for l in poly_name_in_collapse] \
        if poly_name_in_collapse != [] else []

    # Add gap, N, stop codon and mt2var markers
    trim_args = parse_trim_arg(trim_kw)
    add_markers(aln_d, site_type, poly_sample_prefix_list, 
                mt2_seg_allele_filter_mode, only_conserved_aa, genetic_code_type,
                trim_args)

    # Make site position table
    df = get_orig_pos_table(aln_d, site_type)

    # Output unfiltered alignments and table
    save_to_FASTA_files(
        aln_d, 
        dir_paths[0], 
        sample_order,
        input_items.prefix,
        input_items.suffix, 
        include_markers=True, output_empty_aln=False
    )

    # Filter alignments and apply filtering to table
    print('\nRetain codons with final marker 1')
    retain_sites_in_AlnSet(aln_d)

    # Output filtered alignments and table
    save_to_FASTA_files(
        aln_d, 
        dir_paths[1], # Filtered codon/intron alignment
        sample_order,
        input_items.prefix,
        input_items.suffix, 
        include_markers=False, output_empty_aln=False
    )
    df.to_csv(dir_paths[1]+'.csv', index=False)

    # Extract 3rd positions
    if site_type == 'codon':
        print('\nExtract 3rd positions')
        retain_site_from_codon_in_AlnSet(aln_d, 3)
        save_to_FASTA_files(
            aln_d, 
            dir_paths[2], # Filtered codon 3rd pos alignments
            sample_order,
            input_items.prefix,
            input_items.suffix, 
            include_markers=False, output_empty_aln=False
        )

    # # Make collapse sequences
    # print('\nMake collapse sequneces')
    # # poly_name_in_collapse=[['RG', 'RG_collapse'], ['MD', 'MD_collapse']]
    # if poly_sample_prefix_list != []:
    #     print(poly_info)
    #     make_collapse_sequences(aln_d, poly_name_in_collapse)
    #     collapse_sample_order = [
    #         sample
    #         for sample in sample_order
    #         if len([p for p in poly_info 
    #                 if not sample.startswith(p.sample_prefix)]) == len(poly_info)
    #     ] + [p.collapse_prefix+s for p in poly_info for s in ['_a', '_b']]
    #     print(collapse_sample_order)
    #     # print(aln_d[1].sample_ids)

    #     save_to_FASTA_files(
    #         aln_d, 
    #         dir_paths[-2], 
    #         collapse_sample_order,
    #         input_items.prefix,
    #         input_items.suffix, 
    #         include_markers=False, output_empty_aln=False
    #     )

    # # Concatenate non-collapse sequences
    # dir_to_concat = dir_paths[2] if site_type == 'codon' else dir_paths[1]
    # print('\nConcatenate sequences')
    # concatenate_sequences(
    #     os.path.join(dir_paths[-1], 'concat_non-collapse.aln'), # out path
    #     dir_to_concat, # input aln dir
    #     os.path.join(dir_to_concat, '0.itemlist'), # item list
    #     assert_file_exist=assert_file_exist
    # )
    
# ------- sub-main functions ------- #

def check_output_paths(aln_out_dir, site_type):
    if site_type == 'codon':
        dir_names = [
            '1_unfiltered_with_markers',
            '2_filtered_without_markers', # with table
            '3_filtered_without_markers_3rd_pos', 
            # '4_collapse_seqs_3rd_pos', 
            # '5_non_collapse_concat_seq_3rd_pos'
        ]
    elif site_type == 'nucleotide':
        dir_names = [
            '1_unfiltered_with_markers',
            '2_filtered_without_markers', # with table
            # '3_collapse_seqs', 
            # '4_non_collapse_concat_seq'
        ]

    dir_paths = [os.path.join(aln_out_dir, name) for name in dir_names]
    if not os.path.isdir(aln_out_dir):
        os.mkdir(aln_out_dir)
    for dir_path in dir_paths:
        os.mkdir(dir_path)

    if os.path.isfile(dir_paths[1]+'.csv'):
        raise FileExistsError(dir_paths[1]+'.csv')

    return dir_paths

def check_sample_size(aln_d):
    nsamples = [aln.nsamples for k, aln in aln_d.items()]
    if min(nsamples) != max(nsamples):
        raise ValueError(
            f'Number of samples are different somewhere: {min(nsamples)} != {max(nsamples)}')

def read_alignments_from_directory(
        input_aln_dir, input_aln_list_file, marker_kw, markers_to_use,
        assert_file_exist):

    # Read item list and get file list
    fpath_list, items = get_file_list_from_item_list_file(
        input_aln_list_file, input_aln_dir)
    print(len(fpath_list), 'alignment files are found.')

    # Read alignment file to a dictionary
    filename_to_key_encoder = lambda fpath: os.path.basename(fpath)\
        .split(items.prefix)[1]\
        .split(items.suffix)[0]
    aln_d = read_fasta_files_to_Alignment_d(
        fpath_list, marker_kw, markers_to_use,
        assert_file_exist, filename_to_key_encoder)
    print(len(aln_d), 'alignments are read.')
    return aln_d, items

def add_markers(
        aln_d:dict, # A dictionary containing Alignment objects
        site_type:str, # 'codon' or 'nucleotide'
        poly_sample_prefix_list:list, # A list of population sample prefix names
        mt2_seg_allele_filter_mode:str, # 'any' or 'third_pos' (CDS alignment only)
        only_conserved_aa:bool, # Whether to filter non-conserved amino acid (CDS alignment only)
        genetic_code_type:str, # Type of GENETIC CODE. "a" or "b"
        trim_args:dict={},
        ):
    """ Add markers of gap, N, stop codon and more than two variants across
    codons. All sample sequences are used. """
    if genetic_code_type == 'a':
        genetic_code_d = GENETIC_CODE_a
    elif genetic_code_type == 'b':
        genetic_code_d = GENETIC_CODE_b
    else:
        raise TypeError(f'Unknown genetic code type "{genetic_code_type}" is '\
                         'found. Please specify "a" or "b".')

    print('\nAdd markers')
    # An option whether to filter codons for non-conserved amino acids is 
    # applied only for codon alignment but not for nucleotide alignments like 
    # intron sequence alignment.
    only_conserved_aa = False if site_type == 'nucleotide' else only_conserved_aa
    for cnt, item in enumerate(aln_d.items()):
        _, aln = item
        # Add gap marker: nucleotide (intron) is supported.
        _add_marker_to_Aln(
            aln, filt_char='-', marker_name='gap_marker', site_type=site_type)

        # Add N marker: nucleotide (intron) is supported.
        _add_marker_to_Aln(
            aln, filt_char='N', marker_name='N_marker', site_type=site_type)

        # Add stop codon marker
        if site_type == 'codon':
            _add_stopcodon_marker_to_Aln(aln)
            
        # Add multistate marker
        if poly_sample_prefix_list != []:
            _add_mt2var_marker_to_Aln(
                aln=aln, 
                site_type=site_type, # nucleotide (intron) is supported.
                poly_sample_prefix_list=poly_sample_prefix_list,
                mode=mt2_seg_allele_filter_mode
            )

        # Add trim marker if arguments are given
        if len(trim_args) > 0:
            _add_trim_marker(aln, site_type, **trim_args)
        
        if only_conserved_aa:
            _add_consAA_marker_to_Aln(aln, genetic_code_d=genetic_code_d)
            # nucleotide (intron) is supported.
            _add_final_marker(aln, site_type, aa_marker_name='AA_marker')

        else:
            # nucleotide (intron) is supported.
            _add_final_marker(aln, site_type, aa_marker_name='')

        if (cnt+1) % 1000 == 0:
            print(cnt+1)
            
    print(cnt+1)

def retain_site_from_codon_in_AlnSet(aln_d, position):
    for cnt, item in enumerate(aln_d.items()):
        _, aln = item
        retain_site_from_codon_in_Aln(aln, position)
        
        if (cnt+1) % 1000 == 0:
            print(cnt+1)
            
    print(cnt+1)

def make_collapse_sequences(
    aln_d, 
    poly_name_in_collapse, 
    suffix=('a', 'b')):
    for cnt, item in enumerate(aln_d.items()):
        geneinfo_id, aln = item
        if aln.nsites == 0:
            continue
        
        make_collapse_aln(aln, poly_name_in_collapse, suffix)

        if cnt % 1000 == 0:
            print(cnt)
    print(cnt)

def concatenate_sequences(
    out_path, input_aln_dir, input_aln_list_file, assert_file_exist):
    aln_d, items = read_alignments_from_directory(
        input_aln_dir, input_aln_list_file, marker_kw=None, markers_to_use=[],
        assert_file_exist=assert_file_exist)
    
    # Concatenate alignments
    concat = concatenate_alignments(aln_d, items.item_list)
    fasta_lines = [
        f'>{name}\n{seq}'
        for name, seq in concat.items()
    ]

    with open(out_path, 'w') as f:
        print('\n'.join(fasta_lines), file=f)

    to_item_list(
        os.path.join(os.path.dirname(out_path), 'baseml_bin_list'), # path
        'concat_collapse', # name
        items.item_list, # list
        items.prefix, # prefix
        items.suffix # suffix
    )

# ------- small functions ------- #

def concatenate_alignments(aln_d, order):
    concat_names = [a.split('_')[0] for a in aln_d[order[0]].sample_ids]
    concat_seqs = np.array([
        aln_d[name].sample_sequences
        for name in order
    ])

    concat_seqs_t = [
        ''.join(seqs)
        for seqs in concat_seqs.T
    ]

    return dict(zip(concat_names, concat_seqs_t))

def collapse(site_str):
    unique_sites = list(set(site_str))
    if len(unique_sites) == 1:
        return [unique_sites[0], unique_sites[0]]
    shuffle(unique_sites)
    return unique_sites

def make_collapse_aln(
    aln, poly_name_in_collapse, suffix=('a', 'b')):
    # Collapse sequences for Dmel and Dsim alleles
    samples_to_remove = []

    for key, out_prefix in poly_name_in_collapse:
        # Generate a new alignment containing only samples with names
        # starting with Dmel_RG
        aln_subset = aln.get_samples([key], match_prefix=True).samples
        samples_to_remove += aln_subset.ids

        # Create a numpy array of the sequences and transpose it such that
        # the sites become rows
        # NOTE: Do not use transpose because it will take forever.
        aln_subset_tmat = np.array([list(s) for s in aln_subset.sequences]).T

        # Collapse each row into two columns by 
        aln_collapsed = [collapse(site) for site in aln_subset_tmat]
        aln_coll_mat = [''.join(row) for row in np.array(aln_collapsed).T]
        aln.append_sample_from_lists(
            ['{}_{}'.format(out_prefix, suffix[0]), 
             '{}_{}'.format(out_prefix, suffix[1])],
            ['', ''],
            aln_coll_mat
        )

    # Remove Dmel_RG and Dsim_MD alleles
    aln.remove_samples(samples_to_remove, match_prefix=False)

def retain_site_from_codon_in_Aln(aln, position):
    assert position in {1, 2, 3}, '1, 2 or 3 is acceptable for position.'
    slide = position-1
    keep_pos_list = [
        cod+slide for cod in range(0, aln.nsites, 3)
    ]
    
    aln.retain_sites(keep_pos_list)

def get_orig_pos_table(aln_d, site_type):
    cols = ['aln_name', 'orig_site_pos']
    concat_list = [
        get_pos_table(aln).assign(aln_name=name)
        for name, aln in aln_d.items()
    ]
    if site_type == 'codon':
        df = pd\
            .concat(concat_list)\
            .reset_index()\
            .rename(columns={'index': 'orig_site_pos'})\
            .assign(thirdpos=lambda d: d['orig_site_pos']\
                .apply(lambda s: 1 if s % 3 == 2 else 0))

        return df[(df['final_marker'] == 1) & (df['thirdpos'] == 1)]\
            .loc[:, cols]
        
    elif site_type == 'nucleotide':
        df = pd\
            .concat(concat_list)\
            .reset_index()\
            .rename(columns={'index': 'orig_site_pos'})

        return df[(df['final_marker'] == 1)].loc[:, cols]
        
    raise TypeError(f'Unexpected site_type {site_type} is found. '\
                     'Please specify "codon" or "nucleotide".')

def get_pos_table(aln):
    marker_name = 'final_marker'
    marker_seqs = aln.get_markers([marker_name]).sequences
    assert len(marker_seqs) == 1
    marker_seq = marker_seqs[0]
    return pd.DataFrame({marker_name: list(map(int, marker_seq))})

def add_gap_marker_to_AlnSet(aln_d, site_type):
    # Add gap marker
    for i, item in enumerate(aln_d.items()):
        _, aln = item
        
        _add_marker_to_Aln(
            aln, filt_char='-', marker_name='gap_marker', site_type=site_type)
        
        if (i+1) % 1000 == 0:
            print(i+1)
    print(i+1)

def add_N_marker_to_AlnSet(aln_d, site_type):
    # Add gap marker
    for i, item in enumerate(aln_d.items()):
        _, aln = item
        
        _add_marker_to_Aln(
            aln, filt_char='N', marker_name='N_marker', site_type=site_type)
        
        if (i+1) % 1000 == 0:
            print(i+1)
    print(i+1)

def _add_marker_to_Aln(aln, site_type, filt_char, marker_name):
    size = get_column_size(site_type)
    
    # all samples registered in the aln will be used
    subset_samples = aln.sample_ids

    # Create an initial filter array of 1
    # Length is number of codons or sites
    filter_array = np.ones(int(aln.nsites/size))
    # Determine sites with N in any codon
    position_list = [
    i
    # Loop over sample sites by 3's, codon_sites is a list of 3-char strings
    for i, sites in enumerate(aln.get_samples(subset_samples)\
                                 .iter_sample_sites(size=size))
        # Loop over each unique variant of 3-char strings
        for variant in set(sites)
        # If "N" character is found, include the current position i
        if filt_char in variant.upper()
    ]
    filter_array[position_list] = 0
    
    # Add new marker
    aln.markers.append_rows(
        [marker_name],
        ['notes="{} if codon site has no {} char, else {}"'\
            .format('1'*size, filt_char, '0'*size)],
        [''.join(['1'*size if i else '0'*size for i in filter_array])]
    )

def add_stopcodon_marker_to_AlnSet(aln_d):
    # Add stop codon marker
    for i, item in enumerate(aln_d.items()):
        _, aln = item
        _add_stopcodon_marker_to_Aln(aln)

        if (i+1) % 1000 == 0:
            print(i+1)
    print(i+1)

def _add_stopcodon_marker_to_Aln(aln):
    """ Maybe not necessary because user will specify 
    individual amino acid for final alignments """
    # Create an initial filter array of 1
    # Length is number of codons
    filter_array = np.ones(int(aln.nsites/3))
    subset_samples = aln.sample_ids

    # Determine sites with gap "-" in any codon
    position_list = [
        i
        # Loop over sample sites by 3's, codon_sites is a list of 3-char strings
        for i, codon_sites in enumerate(aln.get_samples(subset_samples)\
                                           .iter_sample_sites(size=3))
        # Loop over each unique variant of 3-char strings
        for variant in set(codon_sites)
        # If "-" character is found, include the current position i
        if variant.upper() in STOPCODONS
    ]
    filter_array[position_list] = 0
    
    # Add new marker
    aln.markers.append_rows(
        ['stop_codon_marker'],
        ['notes="111 if codon site is not a stop codon, else 000"'],
        [''.join(['111' if i else '000' for i in filter_array])]
    )

def add_mt2var_marker_to_AlnSet(aln_d, 
    site_type, # 'codon' or 'nucleotide' is accepted
    poly_sample_prefix_list, # a list of sample name prefix to indicate population samples
    mode='third_pos',
    ):
    # Add trichar at any site marker
    for i, item in enumerate(aln_d.items()):
        _, aln = item
        
        _add_mt2var_marker_to_Aln(
            aln=aln, 
            site_type=site_type, 
            poly_sample_prefix_list=poly_sample_prefix_list,
            mode=mode
        )
        
        if (i+1) % 1000 == 0:
            print(i+1)
    print(i+1)

def _add_mt2var_marker_to_Aln(
    aln, # Alignment object
    site_type:str, # 'codon' or 'nucleotide' is accepted
    poly_sample_prefix_list, # a list of sample name prefix to indicate population samples
    mode='third_pos', # If "any" is given, then a codon with more than three 
    # vars at any site will get 0
    ):
    
    size = get_column_size(site_type)
    # Overwrite mode argument if site_type is nucleotides
    mode = 'nucleotide' if site_type == 'nucleotide' else mode
    
    # These characters will be ignored from counting of variants
    avoid_chars = {'N', '-'} 
    multiple_var_sites_in_codon = lambda sites_set: \
        len([
            n for n in range(3) 
            if len(set([codon[n] for codon in sites_set])) >= 2]) \
        >= 2 # True if multiple sites are segregating in a pop
    
    for n, sample_id_prefix in enumerate(poly_sample_prefix_list):
        # Initialize codon/site array
        # [1, 1, 1] represents positions of three codons but not three 
        # nucleotides from a single codon, when site_type is 'codon'
        filter_array = np.ones(int(aln.nsites/size))

        # A codon/site column get 0, when any of tri nucleotides has more than 
        # two variants within a specified population.
        # filter_pos_list will contains positions of codon/site columns which 
        # will be removed.
        sites_set_generator = (
            set(sites) # set of aligned codons in a population
            for sites in aln.get_samples(sample_id_prefix, match_prefix=True)\
                            .iter_sample_sites(size=size)
        )
        if mode == 'third_pos':
            filter_pos_list = list(set([
                i
                for i, sites_set in enumerate(sites_set_generator)
                if len(set([sites[2] for sites in sites_set 
                    if sites[2] not in avoid_chars])) > 2
            ]))
        elif mode == 'codon':
            # sites_set like {'GGA', 'GGG', 'GGT'} and {'AAA', 'ATT'} are filtered.
            filter_pos_list = list(set([
                i
                for i, sites_set in enumerate(sites_set_generator)
                if (len(sites_set) > 2) or (multiple_var_sites_in_codon(sites_set))
            ]))

        elif mode == 'nucleotide':
            if site_type == 'codon':
                raise TypeError('Please specify "third_pos" or "codon" as mode '\
                                'for codon alignments.')
            # sites_set like {"A", "T", "C"} are filtered
            filter_pos_list = list(set([
                i
                for i, sites_set in enumerate(sites_set_generator)
                if (len(sites_set) > 2)
            ]))

        else:
            raise TypeError(f'Unknown mode is found: {mode}. Please specify '\
                             '"third_pos", "codon" or "nucleotide".')

        # Replace 1 with 0 at a codon/site columns which will be removed.
        filter_array[filter_pos_list] = 0

        # Add new marker
        aln.markers.append_rows(
            ['mt2var_{}_marker'.format(sample_id_prefix)],
            ['notes="{} if sites is not "{}", has 2 or less variants, '\
                'and distance is 1 change away; else {}"'\
                    .format('1'*size, '-'*size, '0'*size)],
            [''.join(['1'*size if i == 1 else '0'*size for i in filter_array])]
        )

def _add_trim_marker(aln, site_type, ref_sample_name, left_num, right_num):
    size = get_column_size(site_type)
    filter_array = np.ones(int(aln.nsites/size))

    # These characters will be ignored from counting of variants
    ref_seq = aln.get_samples(ref_sample_name, match_prefix=True).sample_sequences[0]
    aln_len = len(ref_seq)
    ref_seq_len = aln_len - ref_seq.count('-')

    # If ref pos is greater than or equal to this value, the position is filtered.
    right_pos = ref_seq_len+1 - right_num 

    filter_pos_list = []
    pos_in_ref_seq = 0

    # For each site/codon of alignment
    for i in range(0, aln_len-(size-1), size):
        # Get reference sequence
        ref_sites = ref_seq[i:i+size]

        # Check no gap in ref site
        if '-' in set(ref_seq[i:i+size]):
            continue

        pos_in_ref_seq += 1 # The first site is 1

        # Trim left
        if pos_in_ref_seq <= left_num:
            filter_pos_list.append(i)

        # Trim right
        elif pos_in_ref_seq >= right_pos:
            filter_pos_list.append(i)

    # Replace 1 with 0 at a codon/site columns which will be removed.
    filter_array[filter_pos_list] = 0

    # Add new marker
    aln.markers.append_rows(
        ['trim_by_{}_pos_marker'.format(ref_sample_name)],
        ['notes="{} if sites is not "-" and not within trim range; else {}"'\
                .format('1'*size, '0'*size)],
        [''.join(['1'*size if i == 1 else '0'*size for i in filter_array])]
    )

def parse_trim_arg(trim_kw):
    if trim_kw == '':
        return {}

    parts = trim_kw.split(':')
    left, right = parts[1].split('..')[:2]

    return {
        'ref_sample_name': parts[0].strip(),
        'left_num': int(left.strip()),
        'right_num': int(right.strip())
    }

def get_column_size(site_type):
    if site_type == 'codon':
        return 3
    elif site_type == 'nucleotide':
        return 1

    raise Exception('site_type has to be "codon" or "nucleotide".')    

def add_consAA_marker_to_AlnSet(aln_d, genetic_code_d):
    """ Add conserved or non-conserved marker based on keep_cons flag. """
    # Add consAA marker
    for i, item in enumerate(aln_d.items()):
        _, aln = item
        _add_consAA_marker_to_Aln(aln=aln, genetic_code_d=genetic_code_d)

        if (i+1) % 1000 == 0:
            print(i+1)
    print(i+1)

def _add_consAA_marker_to_Aln(aln, genetic_code_d):
    # if empty list is passed to subset_samples, 
    # all samples registered in the aln will be used
    subset_samples = aln.sample_ids
        
    # compare AA across samples
    is_conserved = lambda x: min(x) == max(x)
    
    # Collect aa sequences in aa_matrix_t where rows correspond to sites
    # array([['M', 'M', 'M', ..., 'M', 'M', 'M'],
    #       ...,
    #       ['s', 's', 's', ..., 's', 's', 's']], dsite_type='<U1')
    aa_matrix_t = np.array([
        [genetic_code_d[codon] 
            if ('N' not in codon) and ('-' not in codon) else '_'
            for codon in set(codons)]
                for codons in aln.get_samples(subset_samples, match_prefix=True)\
                                 .iter_sample_sites(size=3)
    ])
    
    # Create filter track that marks conserved sequences as 1
    filter_array = [is_conserved(row) for row in aa_matrix_t]
    
    # Add new marker
    aln.markers.append_rows(
        ['consAA_marker'],
        ['notes="111 if codon site conserved, else 000"'],
        [''.join(['111' if i else '000' for i in filter_array])]
    )
    
    # Add new marker, amino acid identities
    aln.markers.append_rows(
        ['AA_marker'],
        ['notes="Amino acid letter if conserved, else _"'],
        [''.join([str(aa_matrix_t[i][0])*3 if v else '___'
                for i, v in enumerate(filter_array)])]
    )

def _add_final_marker(aln, site_type, aa_marker_name):
    marker_name = 'final_marker'
    if marker_name in aln.marker_ids:
        raise TypeError('"final_marker" exists.')

    size = get_column_size(site_type)

    # Initialize array with a number of codons of 0
    keep_site_array = np.zeros(aln.nsites)

    # Get a list of markers to be used for filtering
    marker_ids = deepcopy(aln.marker_ids)
    # if avoid_AA is True, AA_marker will not be used. 
    # In the future, this will not be removed 
    # from user option and become default
    if aa_marker_name != '':
        try:
            idx = marker_ids.index(aa_marker_name)
        except ValueError:
            # print('"AA_marker" is not found in {}'.format(geneinfo_id))
            idx = slice(0) # nothing will be deleted
        del marker_ids[idx]
    
    # Get marker alignments and turn into a numpy array    
    marker_matrix = np.array([list(map(int, m)) 
        for m in aln.get_markers(marker_ids).sequences])
    # Sum the values down each column
    # Keep positions where sum is equal to the number of rows 
    # of the marker matrix.
    # These columns have passed all filters.
    # Columns whose sum is less than the number of rows 
    # have failed one or more filters
    keep_pos_list = [
        p
        for i in range(0, marker_matrix.shape[1], size)
            if marker_matrix.T[i:i+size].sum() == len(marker_ids) * size
            for p in range(i, i+size)
    ]

    # This assertion is for CDS alignment.
    assert len(keep_pos_list) % size == 0, \
        f'Positions to keep in the final set is not multiple of '\
        f'{size}: {len(keep_pos_list)}'
    keep_site_array[keep_pos_list] = 1
    
    aln.markers.append_rows(
        [marker_name],
        ['notes="111 if none of markers has 000 otherwise 000"'],
        [''.join([str(int(i)) for i in keep_site_array])]
    )
    
def add_final_marker_to_AlnSet(aln_d, site_type, aa_marker_name):
    for cnt, items in enumerate(aln_d.items()):
        _, aln = items
        _add_final_marker(aln, site_type, aa_marker_name)
        
        if (cnt+1) % 1000 == 0:
            print(cnt+1)
    print(cnt+1)

def retain_sites_in_Aln(aln):
    marker_name='final_marker'
    filt_marker_array = np.array(
        list(map(int, aln.get_markers(marker_name).sequences[0])))
    keep_pos_array = np.where(filt_marker_array == 1)[0]
    
    aln.retain_sites(keep_pos_array)
    
def retain_sites_in_AlnSet(aln_d):
    for cnt, items in enumerate(aln_d.items()):
        name, aln = items
        retain_sites_in_Aln(aln)
        
        if (cnt+1) % 1000 == 0:
            print(cnt+1)
            
    print(cnt+1)

ItemList = namedtuple(
    'ItemList', 
    ['name', 'prefix', 'suffix', 'item_list']
)
def read_item_list(item_list_path:str):
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
                parts = l[1:].split()
                list_name = parts[0]
                itemnum = int(parts[1])
                continue

            else:
                item_list.append(l)
    
    assert len(item_list) == itemnum
    return ItemList(list_name, prefix, suffix, item_list)

def to_item_list(out_path:str, name:str, item_list:list, prefix:str, suffix:str):
    out_lines = [
        '1',
        f'prefix: {prefix}',
        f'suffix: {suffix}',
        f'>{name}\t{len(item_list)}'
    ]

    for item in item_list:
        out_lines.append(str(item))

    with open(out_path, 'w') as f:
        print('\n'.join(out_lines), file=f)

def get_file_list_from_item_list_file(item_list_path, dir_path):
    """ Returns file list and item list. """
    items = read_item_list(item_list_path)
    file_list = []

    for item in items.item_list:
        fname = items.prefix + item + items.suffix
        file_list.append(os.path.join(dir_path, fname))
    
    return file_list, items

def read_fasta_files_to_Alignment_d(
        fpath_list:list,
        marker_kw:str,
        markers_to_use:list,
        assert_file_exist:bool,
        filename_to_key_encoder=None):
    """ Returns a dictionary of Alignment object by reading FASTA file. """
    def func(fpath, aln_name, marker_kw, markers_to_use):
        aln = al.fasta_file_to_alignment(
            fpath, 
            filename_to_key_encoder(fpath),
            marker_kw
        )
        markers_to_remove = [
            m for m in aln.marker_ids if m not in markers_to_use]
        not_found_markers = set(markers_to_use) - set(aln.marker_ids)
        if len(not_found_markers) > 0:
            warn(f'{not_found_markers} are not found in the expected markers.')

        # Convert C/N/G marker to 1/0 marker
        if isinstance(markers_to_use, list):
            if len(markers_to_use) > 0:
                markers_to_convert, converted_seqs = zip(*[
                    (name, ''.join(['1' if a == 'C' else '0' for a in seq ]))
                    for name, seq in zip(
                        markers_to_use, aln.get_markers(markers_to_use).sequences)
                    if len(set(seq) - {'C', 'G', 'N'}) == 0 # If sequnece has only C, G or N
                ])
            else:
                markers_to_convert, converted_seqs = [], []
        else:
            raise TypeError('Wrong type of marker list is found: '\
                            f'{type(markers_to_use)} for {markers_to_use}')
        
        # Remove markers that will not be used and that will be replaced.
        markers_to_remove += list(markers_to_convert)
        if len(markers_to_remove) > 0:
            aln.markers.remove_rows_by_name(markers_to_remove)

        # Add converted markers 
        if len(markers_to_convert) > 0:
            aln.markers.append_rows(
                markers_to_convert,
                ['']*len(markers_to_convert),
                converted_seqs
            )
        return aln

    if filename_to_key_encoder == None:
        filename_to_key_encoder = lambda s: s

    if assert_file_exist:
        return {
            filename_to_key_encoder(fpath): func(
                fpath, filename_to_key_encoder(fpath), 
                marker_kw, markers_to_use)
            for fpath in fpath_list
        }
    return {
        filename_to_key_encoder(fpath): func(
            fpath, filename_to_key_encoder(fpath), 
            marker_kw, markers_to_use)
        for fpath in fpath_list
        if os.path.isfile(fpath)
    }

def save_to_FASTA_files(
        aln_d:dict, 
        aln_out_dir:str, 
        sample_order:list,
        prefix:str,
        suffix:str,
        include_markers:bool,
        output_empty_aln:bool,
        output_not_matched:bool=False,
        ):

    saved_count = 0
    item_list = []
    for name, aln in aln_d.items():
        if aln.nsites == 0:
            if not output_empty_aln:
                continue
                
        fasta_path = os.path.join(aln_out_dir, prefix+name+suffix)

        fasta_line = aln_to_FASTA(
            aln, sample_order, include_markers, output_not_matched)
        with open(fasta_path, 'w') as f:
            print(fasta_line, file=f)

        item_list.append(name)
        saved_count += 1

    out_list_path = os.path.join(aln_out_dir, '0.itemlist')
    to_item_list(
        out_list_path, 
        os.path.basename(aln_out_dir.rstrip('/')), 
        item_list, 
        prefix, 
        suffix
    )

def aln_to_FASTA(aln, sample_order, include_markers, output_not_matched=False):
    def sort_sample_ids(aln, sample_order, output_not_matched=False):
        ordered_samples = []
        for sample_query in sample_order:
            ordered_samples += [sample 
                for sample in aln.sample_ids 
                if re.match(sample_query, sample)
            ]

        if output_not_matched:
            ordered_samples += [
                sample
                for sample in aln.sample_ids
                if sample not in ordered_samples
            ]

        return ordered_samples

    ordered_samples = sort_sample_ids(aln, sample_order, output_not_matched)

    fasta_lines = [
        '>{seqname}\n{seq}'.format(
            seqname=sample,
            seq=aln.get_samples([sample], match_prefix=False).sample_sequences[0]
        )
        for sample in ordered_samples
    ]
    if include_markers:
        fasta_lines += [
            '>{seqname}\n{seq}'.format(
                seqname=marker,
                seq=aln.get_markers(marker).sequences[0]
            )
            for marker in aln.marker_ids
        ]

    return '\n'.join(fasta_lines)

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

    # Check sections. Only 'REQUIRED' and 'OPTIONAL' sections will be used.
    assert 'ALIGNMENT' in config.sections(), \
        f'ALIGNMENT section is not found in {control_file}.'
    expected_sections = ['ALIGNMENT']
    not_expected_sections = [
        s for s in config.sections() if s not in expected_sections]
    if len(not_expected_sections) > 1:
        msg = f'Unexpected sections, {", ".join(not_expected_sections)}, '\
              'were found. These are not used in the analysis. '
        warn(msg)

    converters_d = {
        'poly_info': poly_info_parser_from_ctl,
        'only_conserved_aa': lambda s: s == 'True',
        'assert_file_exist': lambda s: s == 'True',
        'sample_order': lambda s: [a.strip() for a in s.split(',')],
        'markers_to_use': lambda s: [] if s == 'None' else [
            a.strip() for a in s.split(',')],
        'trim_kw': lambda s: '' if s == 'None' else s
    }

    flattened  = [
        (opt, converters_d[opt](v)) 
        if opt in converters_d.keys() else (opt, v) 
            for opt, v in config['ALIGNMENT'].items()
    ]

    return dict(flattened)

PopulationSampleInfo = namedtuple(
    'PopulationSampleInfo', [
        'sample_prefix', # Prefix in sample seq
        'collapse_prefix', # Prefix in collapse seq
        'extant_nodes', # Extant node names
        'ancestor_node', # Ancestor node name
        'sample_seq_start', # Position of the first allele in sample alignment.
        'allele_count' # Number of alleles
    ]
)
def poly_info_parser_from_ctl(poly_info_str):
    populations = [
        l.strip() for l in poly_info_str.split('\n') if l.strip() != '']

    poly_info_list = []
    for pop in populations:
        parts = pop.split()
        try:
            poly_info_list.append(
                PopulationSampleInfo(
                    parts[0], 
                    parts[1],
                    [int(a.strip()) for a in parts[2].split(',')],
                    int(parts[3]), 
                    int(parts[4]), 
                    int(parts[5])
                )
            )
        except IndexError:
            raise IndexError(parts)

    return poly_info_list

def to_control_file(output_file, start_date, **kwargs):
    lines = [f'{k} = {v}' for k, v in kwargs.items()]

    with open(output_file, 'w') as f:
        print(f'Run date: {start_date}.', file=f)
        # print(f'Run as a part of CCE pipeline ({__version__})\n', file=f)
        print('\nThe listed parameters are used for the run.\n', file=f)
        print('\n'.join(lines), file=f)

if __name__ == '__main__':
    control_file = sys.argv[1]
    params = read_control_file(control_file)
    print('\nINPUT')
    for k, v in params.items():
        print(f'{k}: {v}')
    
    main(
        params['input_aln_dir'],
        params['input_aln_list_file'],
        params['aln_out_dir'],
        params['site_type'],
        params['sample_order'],
        params['marker_kw'],
        params['markers_to_use'],
        params['poly_info'],
        params['mt2_seg_allele_filter_mode'],
        params['only_conserved_aa'],
        params['genetic_code_type'],
        params['assert_file_exist'], 
        params['trim_kw']
    )
