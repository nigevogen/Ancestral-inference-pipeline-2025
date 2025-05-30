
import os
import argparse
import configparser
from collections import namedtuple

def main(input_ctl_file, output_ctl_file):
    if os.path.isfile(output_ctl_file):
        raise FileExistsError(output_ctl_file)
    
    args_d = read_control_file(input_ctl_file)
    output_to_shell_control_file(args_d, output_ctl_file)

def output_to_shell_control_file(args_d, output_ctl_file):
    output_args = [
        'output_dir', 
        'sample_tree_file', 
        'site_type',
        'iteration_num', 
        'internal_node_num',
        'mutation_category_num'
    ]
    output_lines = [
        # Rename to ppda_output_dir to avoid overlap with other output_dir 
        # variable in other codes.
        f'ppda_output_dir={args_d[arg]}'
            if arg == 'output_dir' else f'{arg}={args_d[arg]}' 
        for arg in output_args
    ]
    # Parse polymorphism info
    poly_info_list = poly_info_parser_from_ctl(args_d['poly_info'])

    # If any polymorphism is given,
    if len(poly_info_list) > 0:
        poly_args = ['poly_sp_num={}'.format(len(poly_info_list))]

        # Add to polymorphism info arguments to output arguments
        for n, poly_info in enumerate(poly_info_list):
            poly_args.append(f'poly_sp{n+1}_sample_start={poly_info[4]+1}')
            poly_args.append(f'poly_sp{n+1}_sample_end={poly_info[4]+poly_info[5]}')
            poly_args.append(f'poly_sp{n+1}_collapse1={poly_info[2][0]}')
            poly_args.append(f'poly_sp{n+1}_collapse2={poly_info[2][1]}')
    else:
        raise Exception('No polymorphism info is given. Please check '\
                       f'poly_info: {args_d["poly_info"]}.')

    output_lines += poly_args
    
    # Parse single allele species
    single_allele_sp_list = single_allele_sp_parse_from_ctl(args_d['single_allele_sp'])
    sa_args = ['single_allele_sp_num={}'.format(len(single_allele_sp_list))]

    # Add to single allele arguments to output arguments
    for m, single_allele in enumerate(single_allele_sp_list):
        sa_args.append(f'single_allele_sp{m+1}_sample_pos={single_allele[0]+1}')
    
    output_lines += sa_args
    output_lines.append('total_seq_num_in_collapse={}'\
        .format(len(poly_info_list)*2+len(single_allele_sp_list)))

    with open(output_ctl_file, 'w') as f:
        print('\n'.join(output_lines), file=f)

def poly_info_parser_from_ctl(poly_info_str):
    populations = [
        l.strip() for l in poly_info_str.split('\n') if l.strip() != '']

    poly_info_list = []
    for pop in populations:
        parts = pop.split()
        try:
            poly_info_list.append(
                (
                    parts[0], # Prefix in sample seq
                    parts[1], # Prefix in collapse seq
                    [int(a.strip()) for a in parts[2].split(',')], # Extant node names
                    int(parts[3]), # Ancestor node name
                    int(parts[4]), # Position of the first allele in sample alignment.
                    int(parts[5]) # Number of alleles
                )
            )
        except IndexError:
            raise IndexError(parts)

    return poly_info_list

def single_allele_sp_parse_from_ctl(single_allele_sp_str):
    parts = [a.strip() for a in single_allele_sp_str.split(',')]
    if len(parts) == 0:
        return []

    return [
        tuple([int(a) for a in part.split(':')])
        for part in parts
    ]

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
        arg: value
        for section in expected_sections
            for arg, value in config[section].items()
    }


if __name__ == '__main__':
    desc = 'Make a control file, from which shell variables can be loaded '\
           '`by . <filename>`.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "input_control_file", help="A path to an input control file for pipeline.")
    parser.add_argument(
        "output_control_file", help="A path to an output control file.")

    # Parse arguments from terminal
    args = parser.parse_args()
    main(
        args.input_control_file, 
        args.output_control_file 
    )
