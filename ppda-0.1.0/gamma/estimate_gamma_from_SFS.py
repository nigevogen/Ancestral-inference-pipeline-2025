
import os
import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from popgen.classes import SFS
from popgen.formatter import rewrite_glemin_ctl_template
from popgen.formatter import get_Glemin_maxlnL_bootstrap_results_df

def main(
        gamma_out_dir, 
        test_sfs_dir, 
        test_WS_class_mutation,
        test_SW_class_mutation,
        sfs_anc_node_name, 
        neutral_sfs_dir, 
        sample_allele_num,
        ctl_template_file,
        constraint,
        m=1, 
        theta_range='from_neutral_sfs', 
        r_range=(0.001, 100),
        neutral_mutations=['TA', 'AT', 'GC', 'CG'],
        rep_num=0, 
        rep_num_start=0, 
        rep_num_end=0,
        output_maxlnL_table=True
        ):
    start = datetime.now().isoformat()

    # Make output directories
    make_out_dirs(gamma_out_dir)

    # Check if the final output file does not exist yet.
    maxlnL_table_file = os.path.join(gamma_out_dir, 'maxlnL_result.csv')
    if output_maxlnL_table:
        if os.path.isfile(maxlnL_table_file):
            raise FileExistsError(maxlnL_table_file)

    # Check if constraint argument is valid.
    assert constraint in {'none', 'M1'}, \
        f'Unknown constraint {constraint} is found. Please input "none" or "M1".'

    if rep_num_end == 0:
        rep_num_end = rep_num_start+rep_num

    # Output read arguments as a text file
    run_info_file = os.path.join(gamma_out_dir, '0.run_info.txt')
    to_control_file(
        run_info_file, start,
        gamma_out_dir=gamma_out_dir, 
        test_sfs_dir=test_sfs_dir, 
        test_WS_class_mutation=test_WS_class_mutation,
        test_SW_class_mutation=test_SW_class_mutation,
        sfs_anc_node_name=sfs_anc_node_name, 
        neutral_sfs_dir=neutral_sfs_dir, 
        sample_allele_num=sample_allele_num,
        rep_num=rep_num, 
        ctl_template_file=ctl_template_file,
        constraint=constraint,
        m=m, 
        theta_range=theta_range, 
        r_range=r_range,
        neutral_mutations=neutral_mutations,
        rep_num_start=rep_num_start,
        rep_num_end=rep_num_end
    )

    # Get SFS data of test and neutral classes for all replicates
    if rep_num_end != 0:
        rep_num = rep_num_end - rep_num_start + 1
    else:
        if rep_num == 0:
            raise ValueError('Please input rep_num or rep_num_end.')
        # To run rep_num 0<= rep <= 300, set rep_num_start = 0 and rep_num = 301
        rep_num_end = rep_num_start + rep_num - 1

    print(f'Read {test_sfs_dir}')
    test_sfs = get_sfs_obj(
        test_sfs_dir, sfs_anc_node_name, sample_allele_num, 
        rep_num_start, rep_num_end)
    print(f'Read {neutral_sfs_dir}')
    neutral_sfs = get_sfs_obj(
        neutral_sfs_dir, sfs_anc_node_name, sample_allele_num, 
        rep_num_start, rep_num_end)

    # Make control file and estimate gamma for all replicates
    estimate_gamma_for_all_replicates(
        gamma_out_dir, 
        test_sfs, 
        test_WS_class_mutation,
        test_SW_class_mutation,
        neutral_sfs, 
        int(rep_num_start), 
        int(rep_num_end), 
        ctl_template_file,
        constraint,
        m, 
        theta_range, 
        r_range,
        neutral_mutations
    )

    # Collect estimated gamma with maximum likelihood and make table
    if output_maxlnL_table:
        out_df = get_Glemin_maxlnL_bootstrap_results_df(
            os.path.join(gamma_out_dir, 'rep_out'))

        # Save as a CSV file
        out_df.to_csv(maxlnL_table_file, index=False)

def make_out_dirs(gamma_out_dir):
    if not os.path.isdir(gamma_out_dir):
        os.mkdir(gamma_out_dir)
    if not os.path.isdir(os.path.join(gamma_out_dir, 'rep_ctl')):
        os.mkdir(os.path.join(gamma_out_dir, 'rep_ctl'))
    if not os.path.isdir(os.path.join(gamma_out_dir, 'rep_log')):
        os.mkdir(os.path.join(gamma_out_dir, 'rep_log'))
    if not os.path.isdir(os.path.join(gamma_out_dir, 'rep_out')):
        os.mkdir(os.path.join(gamma_out_dir, 'rep_out'))

def get_sfs_obj(
        sfs_dir, 
        sfs_anc_node_name, 
        sample_allele_num, 
        rep_num_start, 
        rep_num_end):
    concat_list = []
    cols = [
        'anc_node', 'der_node', 'anc_state', 'der_state', 
        'rep', 'frequency', 'count'
    ]
    for n in range(rep_num_start, rep_num_end+1):
        tmp_df = pd\
            .read_csv(os.path.join(sfs_dir, f'rep_{n}/mutation_count.csv'))\
            .assign(rep=n)
        concat_list.append(
            tmp_df[(tmp_df['anc_node'] == sfs_anc_node_name) | 
                   (tmp_df['der_node'] == sfs_anc_node_name)].loc[:, cols]
        ) 
        
    group_cols = ['anc_node', 'anc_state', 'der_state', 'rep', 'frequency']
    df = pd.concat(concat_list)\
        .groupby(group_cols)['count']\
        .sum()\
        .to_frame()\
        .reset_index()

    return SFS(df, sample_allele_num)

def estimate_gamma_for_all_replicates(
        gamma_out_dir, 
        test_sfs, 
        test_WS_class_mutation,
        test_SW_class_mutation,
        neutral_sfs, 
        rep_num_start, 
        rep_num_end, 
        ctl_template_file,
        constraint,
        m=1, 
        theta_range='from_neutral_sfs', 
        r_range=(0.001, 100),
        neutral_mutations=['TA', 'AT', 'GC', 'CG']
        ):
    assert test_sfs.allele_num == neutral_sfs.allele_num, \
        'Wrong allele_num between test and neutral: '\
        f'{test_sfs.allele_num} != {neutral_sfs.allele_num}.'
    allele_num = test_sfs.allele_num

    neu_mutation_t = [tuple(m) for m in neutral_mutations]
    ws_mutation_pool = ['TC', 'TG', 'AG', 'AC']
    sw_mutation_pool = ['CT', 'GT', 'GA', 'CA']
    ws_mutation_pool_2f = ['TC', 'AG']
    sw_mutation_pool_2f = ['CT', 'GA']
    # Set WS mutations
    if test_WS_class_mutation == 'WS':
        ws_mut = [tuple(m) for m in ws_mutation_pool]
    elif test_WS_class_mutation == 'WS2f':
        ws_mut = [tuple(m) for m in ws_mutation_pool_2f]
    else:
        ws_mut = tuple(test_WS_class_mutation)
    # Set SW mutations
    if test_SW_class_mutation == 'SW':
        sw_mut = [tuple(m) for m in sw_mutation_pool]
    elif test_SW_class_mutation == 'SW2f':
        sw_mut = [tuple(m) for m in sw_mutation_pool_2f]
    else:
        sw_mut = tuple(test_SW_class_mutation)

    for n in range(rep_num_start, rep_num_end+1):
        print(n)
        out_ctl_path = os.path.join(gamma_out_dir, f'rep_ctl/ctl_rep_{n}.txt')
        
        neu_table = neutral_sfs.get_table(
            filter_kw={'rep': n, 'frequency': f'ne{allele_num}'},
            groupby_kw={'by': ['anc_state', 'der_state', 'frequency']}, 
            func=np.sum
        ).fillna(0)
        neu = neu_table.loc[neu_mutation_t].sum()
        
        test_table = test_sfs.get_table(
            filter_kw={'rep': n, 'frequency': f'ne{allele_num}'},
            groupby_kw={'by': ['anc_state', 'der_state', 'frequency']}, 
            func=np.sum
        ).fillna(0)
        if test_WS_class_mutation == 'WS' or test_WS_class_mutation == 'WS2f':
            ws = test_table.loc[ws_mut].sum()
        else:
            ws = test_table.loc[ws_mut]
        if test_SW_class_mutation == 'SW' or test_SW_class_mutation == 'SW2f':
            sw = test_table.loc[sw_mut].sum()
        else:
            sw = test_table.loc[sw_mut]
        
        rewrite_glemin_ctl_template(
            pd.DataFrame({'WS': ws, 'SW': sw, 'neu': neu}), 
            ctl_template_file, 
            out_ctl_path, 
            allele_num,
            constraint,
            m, 
            theta_range, 
            r_range
        )

        command = f'anavar '\
            f'"{gamma_out_dir}/rep_ctl/ctl_rep_{n}.txt" '\
            f'"{gamma_out_dir}/rep_out/out_rep_{n}.txt" '\
            f'"{gamma_out_dir}/rep_log/log_rep_{n}.txt"'

        os.system(command)
        
def to_control_file(output_file, start_date, **kwargs):
    lines = [f'{k} = {v}' for k, v in kwargs.items()]

    with open(output_file, 'w') as f:
        print(f'Process start date: {start_date}.', file=f)
        print('Created by estimate_gamma_from_SFS.py\n', file=f)
        print('The listed parameters were used for the run.\n', file=f)
        print('\n'.join(lines), file=f)

if __name__ == '__main__':
    desc = 'Estimate gamma (Ns or selection intensity) from SFS for two '\
           'mutation classes by comparing agains neutral class. This employes '\
           '`anavar` command for maximum likelihood estimation. So `anavar` '\
           'has to be installed somewhere under PATH.'

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        "gamma_out_dir", 
        help='A path to output directory. If a given '\
             'directory does not exist, a new directory will be made.'
    )
    parser.add_argument(
        'test_sfs_dir', 
        help='A path to a directory including SFS for test classes across '\
             'replicates.'
    )
    parser.add_argument(
        'test_WS_class_mutation', help='WS mutation class for test SFS.'
    )
    parser.add_argument(
        'test_SW_class_mutation', help='SW mutation class for test SFS.'
    )
    parser.add_argument(
        'sfs_anc_node_name', 
        help='Ancestor node name for a population of test and neutral SFS',
        type=int
    )
    parser.add_argument(
        'neutral_sfs_dir', 
        help='A path to a directory including SFS for neutral classes across '\
             'replicates.'
    )
    parser.add_argument(
        'sample_allele_num', 
        help='Total number of alleles examined for SFS from the population.',
        type=int
    )
    parser.add_argument(
        "ctl_template_file", 
        help='A path to a template control file for gamma estimation program, '\
             '`anavar`.'
    )
    parser.add_argument(
        'constraint', 
        help='Model constraint between WS and SW mutation classes. This has '\
             'to be "M1" or "none", which are symmetric and asymmetric models, '\
             'respectively.'
    )
    parser.add_argument(
        'rep_num', 
        help='Number of replicates. This has to be the same for test and '\
             'neutral SFS.',
        type=int, nargs='?', default=0
    )
    parser.add_argument(
        'rep_num_start', 
        help='Starting number of replicates. This has to be the same for test and '\
             'neutral SFS. (default: 0)',
        type=int, nargs='?', default=0
    )
    parser.add_argument(
        'rep_num_end', 
        help='Ending number of replicates. This has to be the same for test and '\
             'neutral SFS. If 0 is given, gamma for all replicates are estimated.'\
             '(default: 0)',
        type=int, nargs='?', default=0
    )
    parser.add_argument(
        'output_maxlnL_table', 
        help='Whether to output maximum likelihood tale after estimation. True '\
             'if 1 is given and False if 0 is given (default: 1).',
        type=int, nargs='?', default=1
    )
    # TODO: Make parser for these arguments and make them accessible from 
    # command line.
    # parser.add_argument(
    #     'm', help='SW mutation class for test SFS.',
    #     nargs='?', default=1
    # )
    # parser.add_argument(
    #     'theta_range', help='SW mutation class for test SFS.',
    #     nargs='?', default='from_neutral_sfs'
    # )
    # parser.add_argument(
    #     'r_range', help='SW mutation class for test SFS.',
    #     nargs='?', default=(0.001, 100)
    # )
    # parser.add_argument(
    #     'neutral_mutations', help='SW mutation class for test SFS.',
    #     nargs='?', default=['TA', 'AT', 'GC', 'CG']
    # )

    # Parse arguments from terminal
    args = parser.parse_args()
    
    main(
        args.gamma_out_dir, 
        args.test_sfs_dir, 
        args.test_WS_class_mutation,
        args.test_SW_class_mutation,
        args.sfs_anc_node_name, 
        args.neutral_sfs_dir, 
        args.sample_allele_num,
        args.ctl_template_file,
        args.constraint,
        m=1, 
        theta_range='from_neutral_sfs', 
        r_range=(0.001, 100),
        neutral_mutations=['TA', 'AT', 'GC', 'CG'], 
        rep_num=args.rep_num, 
        rep_num_start=args.rep_num_start,
        rep_num_end=args.rep_num_end,
        output_maxlnL_table=True if args.output_maxlnL_table == 1 else False
    )
