#!/usr/bin/env python3

""" Make a control file for filtering input alignments by reading a control file
for pipeline. 
"""

import os
import argparse
import configparser
from warnings import warn

def main(input_control_file, output_control_file):
    if os.path.isfile(output_control_file):
        raise FileExistsError(output_control_file)
    
    args_d = read_control_file(input_control_file)
    args_d['poly_info'] = reformat_poly_info(args_d['poly_info'])
    output_as_filtering_control_file(args_d, output_control_file)

def output_as_filtering_control_file(args_d, output_control_file):
    # Get required arguments for alignment filtering
    arguments = [
        'input_aln_dir', 
        'input_aln_list_file',
        'aln_out_dir',
        'site_type',
        'sample_order',
        'marker_kw',
        'markers_to_use',
        'poly_info',
        'mt2_seg_allele_filter_mode',
        'only_conserved_aa',
        'genetic_code_type',
        'assert_file_exist', 
        'trim_kw'
    ]
    output_lines = [
        '{} = {}'.format(
            arg, os.path.join(args_d['output_dir'], '1_alignment_processing'))
        if arg == 'aln_out_dir' else f'{arg} = {args_d[arg]}' 
        for arg in arguments
    ]

    with open(output_control_file, 'w') as f:
        print('[ALIGNMENT]', file=f)
        print('\n'.join(output_lines), file=f)

def reformat_poly_info(poly_info_str):
    poly_info_list = [
        '\t'+l.strip()
        for l in poly_info_str.split('\n')
        if l.strip() != ''
    ]
    return '\n' + '\n'.join(poly_info_list)

def read_control_file(input_control_file):
    # Initialize ConfigParser object
    config = configparser.ConfigParser(
        strict=True,
        comment_prefixes=('#', ';'),
        inline_comment_prefixes=('#', ';')
    )

    # Read control file
    paths = config.read(input_control_file)

    # Check number of read control files.
    if len(paths) == 0:
        raise FileNotFoundError(
        f'{input_control_file} is not found.')
    elif len(paths) > 1:
        raise TypeError(f'Iterable {type(input_control_file)} is given as a control '\
            'file. Only one control file is supported.')

    # Check if all the expected sections exist in the control file.
    expected_sections = ['FILES', 'SETTINGS']
    for section in expected_sections:
        assert section in config.sections(), \
            f'{section} section is not found in {input_control_file}.'
    # Flatten section structure
    return {
        arg: value.strip('"')
        for section in expected_sections
            for arg, value in config[section].items()
    }

if __name__ == '__main__':
    desc = 'Make a control file for alignment filtering by reading a control '\
           'file for the entire pipeline.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "input_control_file", help="A path to an input control file for pipeline.")
    parser.add_argument(
        "output_control_file", 
        help="A path to an output control file, which can be used for an input "\
             "of filter_alignments.py")

    # Parse arguments from terminal
    args = parser.parse_args()
    main(
        args.input_control_file, 
        args.output_control_file 
    )
