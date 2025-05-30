import os
import numpy as np
import pandas as pd
import argparse
import configparser
from collections import namedtuple
from datetime import datetime

def main(
        top_out_dir:str, 
        joint_prob_file:str, 
        position_file:str='', 
        position_type:str='', 
        item_list_path:str='', 
        rep_num:int=0 # Number of bootstrap resampling
        ):

    start = datetime.now().isoformat()
    
    # Check output directories
    if not os.path.isdir(top_out_dir):
        os.mkdir(top_out_dir)
    
    bootstrap_out_dir = os.path.join(top_out_dir, 'bootstrap_results')
    if not os.path.isdir(bootstrap_out_dir):
        os.mkdir(bootstrap_out_dir)

    # Output input arguments file
    ctl_file = os.path.join(top_out_dir, '0.read_ctl.txt')
    to_control_file(
        ctl_file, start, 
        top_out_dir=top_out_dir, 
        position_file=position_file, 
        position_type=position_type, 
        item_list_path=item_list_path, 
        joint_prob_file=joint_prob_file, 
        rep_num=rep_num
    )

    # Read position map file
    pos_map = get_pos_map(
        position_file, position_type, item_list_path, top_out_dir)
    if position_file != '':
        print(f'\nPosition map was read from {position_file}')
    
    # Conduct bootstrap resampling per alignment (e.g. CDS or intron)
    if rep_num == 0:
        bootstrap_replicates = []
    else:
        bootstrap_replicates = bootstrap(
            list(pos_map.concat_pos_d.keys()), rep_num)
    
    # Read joint probability file. 
    joint_prob = read_joint_prob_file(joint_prob_file)
    print(f'Joint probability and its frequency was read from {joint_prob_file}')

    count_mutations_by_AWP_method(
        bootstrap_out_dir, 
        joint_prob, 
        pos_map, 
        bootstrap_replicates
    )

# --- sub-functions --- #
def get_pos_map(position_file, position_type, item_list_path, top_out_dir):
    if position_file == '':
        pos_map = []
    else:
        if position_type == 'table':
            items, concat_pos_df = get_concatenated_pos_df(
                position_file, item_list_path)
            
            fname = os.path.basename(position_file)
            pos_map_file = os.path.join(top_out_dir, f'{fname}.posmap.list')
            if os.path.isfile(pos_map_file):
                raise FileExistsError(pos_map_file)

            to_position_map_file(pos_map_file, items, concat_pos_df)

        elif position_type == 'posmap':
            pos_map_file = position_file

        else:
            raise TypeError(f'Unknown position_type {position_type} is found.')

        pos_map = read_position_map_file(pos_map_file)

    return pos_map

def read_joint_prob_file(file_path):
    """ Read a list of joint probabilities of mutations for each lineage with 
    derived allele frequencies and returns a pandas DataFrame. 

    """
    pos_list = []
    node1_list = []
    node2_list = []
    state1_list = []
    state2_list = []
    freq_list = []
    prob_list = []
    changes = []
    pos_in_concat = ''
    
    with open(file_path, 'r') as f:
        for l in f:
            if l.startswith('>'):
                if changes:
                    for n1, n2, m, freq, p in changes:
                        pos_list.append(pos_in_concat)
                        node1_list.append(n1)
                        node2_list.append(n2)
                        state1_list.append(m[0])
                        state2_list.append(m[1])
                        freq_list.append(freq)
                        prob_list.append(p)
                        
                pos_in_concat = int(l[1:].rstrip())
            else:
                changes = parse_probs(l.rstrip())
                
    for n1, n2, m, freq, p in changes:
        pos_list.append(pos_in_concat)
        node1_list.append(n1)
        node2_list.append(n2)
        state1_list.append(m[0])
        state2_list.append(m[1])
        freq_list.append(freq)
        prob_list.append(p)
    
    return pd.DataFrame(
        {
            'pos': pos_list,
            'anc_node': node1_list,
            'der_node': node2_list,
            'anc_state': state1_list,
            'der_state': state2_list,
            'frequency': freq_list,
            'prob': prob_list
        }
    )

def parse_probs(s):
    """ Parse a line of joint probabilities. """
    parts = s.split('\t')
    return [
        tuple(
            [int(node) for node in part.split(':')[0].split('..')] + \
            [part.split(':')[1]] + \
            [int(part.split(':')[2])] + \
            [float(part.split(':')[3])]
        )
        for part in parts
    ]

def count_mutations_by_AWP_method(
        out_dir:str, 
        joint_prob:pd.DataFrame, 
        pos_map:namedtuple=[], 
        bootstrap_replicates:list=[]):
    """ Count the number of mutations using the AWP (Averaging Weighted by 
    Probability) method. 

    Parameters
    ----------
    out_dir: str
        A path to a directory in which mutation table will be output in a 
        sub-directory.
    joint_prob: pandas.DataFrame
        A table of joint probabilities for mutations at each positions. 
    pos_map: PositionMap (default: [])
        A namedtuple object of PositionMap. This can be used to filter 
        mutations by positions.
    bootstrap_replicates: list (default: [])
        A 2D list of alignment items. Each element should include a list of 
        randomly sampled alignment items.
        
    """
    # Make output directory
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    os.mkdir(os.path.join(out_dir, 'rep_0'))

    # Sum probabilities
    group_columns = [
        'anc_node', 'der_node', 'anc_state', 'der_state', 'frequency'
    ]
    # Sum probabilities across all positions
    if len(pos_map) == 0:
        orig_table = joint_prob\
            .groupby(group_columns)['prob']\
            .sum()\
            .to_frame(name='count')\
            .reset_index()

    # Sum probabilities only positions in a given list
    elif isinstance(pos_map, PositionMap):
        positions = [
            p for pos_list in pos_map.concat_pos_d.values() for p in pos_list
        ]
        orig_table = joint_prob[joint_prob['pos'].isin(positions)]\
            .groupby(group_columns)['prob']\
            .sum()\
            .to_frame(name='count')\
            .reset_index()
    else:
        raise TypeError('Invalid object for position map is given.')

    # Output as a file.
    orig_table.to_csv(
        os.path.join(out_dir, 'rep_0/mutation_count.csv'), 
        index=False
    )

    # If bootstrap samples are given,
    if len(bootstrap_replicates) > 0:
        # Concatenated positions in bootstrap replicates were obtained by accessing 
        # from index.
        joint_prob = joint_prob.set_index('pos')
        orig_var_pos_set = set(joint_prob.index)

        # For each bootstrap set,
        for i, rep_id in enumerate(bootstrap_replicates):
            # Make output directory
            os.mkdir(os.path.join(out_dir, f'rep_{i+1}'))
            
            # Get a list of positions in the concatenated sequences for 
            # a given list of alignment items.
            rep_pos = [
                pos 
                for item in rep_id 
                    for pos in pos_map.concat_pos_d[item]
                    # Avoid onvariable sites because they are not in joint prob 
                    # file
                    if pos in orig_var_pos_set 
            ]
            tmp_d = joint_prob.loc[rep_pos]
            
            div_df = tmp_d\
                .groupby(group_columns)['prob']\
                .sum()\
                .to_frame(name='count')\
                .reset_index()
            
            to_pos_list_file(
                rep_id, # aln_name_list
                [pos_map.concat_pos_d[k] for k in rep_id], # concat_pos_list
                [pos_map.orig_pos_d[k] for k in rep_id], # orig_pos_list
                len(rep_pos), # var_site_num
                # out_file
                os.path.join(out_dir, f'rep_{i+1}/used_site_in_concat.list'), 
                False # write_position
            )

            output_file_path = os.path.join(
                out_dir, f'rep_{i+1}/mutation_count.csv')
            div_df.to_csv(output_file_path, index=False)

def to_pos_list_file(aln_name_list, concat_pos_list, orig_pos_list, 
                     var_site_num, out_file, write_position):

    assert len(aln_name_list) == len(concat_pos_list)
    assert len(aln_name_list) == len(orig_pos_list)
    flattened_positions = [
        pos for pos_list in concat_pos_list for pos in pos_list]

    if write_position:
        lines = [
            '>{item}\n{pos_list}'.format(
                item=str(k),
                pos_list='\n'.join([
                    f'{c}\t{o}' for c, o in zip(c_pos, o_pos)
                ])
            ) 
            for k, c_pos, o_pos in zip(aln_name_list, concat_pos_list, orig_pos_list)
        ]
    else:
        lines = [
            f'>{k}' for k in aln_name_list
        ]
    
    with open(out_file, 'w') as f:
        print('aln_num:', len(aln_name_list), file=f)
        print('total_site_num:', len(flattened_positions), file=f)
        print('var_site_num:', var_site_num, file=f)
        print('\n'.join(lines), file=f)

def to_control_file(output_file, start_date, **kwargs):
    lines = [f'{k} = {v}' for k, v in kwargs.items()]

    with open(output_file, 'w') as f:
        print(f'Run date: {start_date}.', file=f)
        # print(f'Run as a part of CCE pipeline ({__version__})\n', file=f)
        print('\nThe listed parameters are used for the run.\n', file=f)
        print('\n'.join(lines), file=f)

# --- Optional functions --- #
def get_concatenated_pos_df(pos_df_path, item_list_path):
    """ Read position map table file into a pandas DataFrame. """
    df = pd.read_csv(pos_df_path)

    items = read_item_list(item_list_path)
    # used_sites = df[(df['final_marker'] == 1)]
    df = df.assign(aln_name=lambda d: d['aln_name'].astype(str))

    # This will be True if no gene is entirely filtered out.
    assert set(df.drop_duplicates(['aln_name'])['aln_name']) \
        == set(items.item_list)

    return items, df\
        .reset_index(drop=True)\
        .set_index('aln_name')\
        .sort_values('orig_site_pos')\
        .loc[items.item_list]\
        .assign(concat_pos=range(len(df.index)))

ItemList = namedtuple(
    'items', 
    ['name', 'prefix', 'suffix', 'item_list']
)

def read_item_list(item_list_path:str):
    """ Read a list of alignment items which was used for concatenation. """
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
                item_list.append(l)
    
    assert len(item_list) == itemnum
    return ItemList(list_name, prefix, suffix, item_list)

def to_position_map_file(output_file, items, concat_pos_df):
    """ Output as an position map list file. """
    out_lines = []

    for item in items.item_list:
        tmp_d = concat_pos_df.loc[item]
        if isinstance(tmp_d, pd.DataFrame):
            concat_pos_list = tmp_d['concat_pos'].tolist()
            orig_pos_list = tmp_d['orig_site_pos'].tolist()
        else:
            concat_pos_list = [tmp_d['concat_pos']]
            orig_pos_list = [tmp_d['orig_site_pos']]

        out_lines.append(
            '>{aln_name}\n{positions}'.format(
                aln_name=item,
                positions= '\n'.join([
                    f'{concat_pos}\t{orig_pos}'
                    for concat_pos, orig_pos in 
                        zip(concat_pos_list, orig_pos_list)
                ])
            )
        )

    assert len(items.item_list) == len(out_lines)

    with open(output_file, 'w') as f:
        print('prefix:', items.prefix, file=f)
        print('suffix:', items.suffix, file=f)
        print(f'itemnum: {len(items.item_list)}', file=f)
        print('\n'.join(out_lines), file=f)

PositionMap = namedtuple(
    'PositionMap',
    ['prefix', 'suffix', 'concat_pos_d', 'orig_pos_d']
)

def read_position_map_file(pos_file):
    """ Read position map list file. """
    item_names = []
    concat_pos_list = []
    temp_concat_pos_list = []
    orig_pos_list = []
    temp_orig_pos_list = []
    item = ''
    
    with open(pos_file, 'r') as f:
        for l in f:
            if l[:-1] == '':
                continue

            if l.startswith('prefix'):
                prefix = l[:-1].split('prefix:')[1].strip()
                continue
                
            if l.startswith('suffix'):
                suffix = l[:-1].split('suffix:')[1].strip()
                continue
                
            if l.startswith('itemnum'):
                exp_itemnum = int(l[:-1].split('itemnum: ')[1])
                continue
                
            if l.startswith('>'):
                if item != '':
                    item_names.append(item)
                    concat_pos_list.append(temp_concat_pos_list)
                    orig_pos_list.append(temp_orig_pos_list)
                    
                item = l[1:-1]
                temp_concat_pos_list = []
                temp_orig_pos_list = []
                
            else:
                pos = l[:-1].split('\t')
                assert len(pos) == 2, pos
                temp_concat_pos_list.append(int(pos[0]))
                temp_orig_pos_list.append(int(pos[1]))
    
    if item != '':
        item_names.append(item)
        concat_pos_list.append(temp_concat_pos_list)
        orig_pos_list.append(temp_orig_pos_list)
        
    concat_d = dict(zip(item_names, concat_pos_list))
    orig_d = dict(zip(item_names, orig_pos_list))
    
    assert len(concat_d) == exp_itemnum, f'{len(concat_d)} != {exp_itemnum}'
    assert len(orig_d) == exp_itemnum, f'{len(orig_d)} != {exp_itemnum}'
    
    return PositionMap(
        prefix, 
        suffix, 
        concat_d,
        orig_d
    )

def bootstrap(data, rep_num, output='matrix'):
    """ Conduct random sampling of items in a given list and returns 
    a matrix of resampled data by default.

    Parameters
    ----------
    data: list
        array-like data
    rep_num: int
        number of replications. bootstrap will be repeated this number of times.
    output: str, 'mean' or 'matrix'
        output info. if mean was specified, only means of bootstrapped data will 
        be outputted.
    """
    data = np.array(data)
    boot_dat = []

    for _ in range(rep_num):
        index = np.random.choice(len(data), len(data), replace=True)
        index.sort()
        if output == 'index':
            boot_dat.append(index)
            continue

        new_data = data[index]
        if output == 'mean':
            boot_dat.append(np.mean(new_data))
        elif output == 'matrix':
            boot_dat.append(np.array(new_data))

    if len(boot_dat) != rep_num:
        raise Exception('bootstrapped data has excess or lack of data.')

    boot_array = np.array(boot_dat)
    # boot_array.sort()
    return boot_array

# --------------#

if __name__ == '__main__':
    desc = 'Output basic statistics for all sites used for ancestral inference.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-o", "--output_dir", 
        help="A path to output directory."
    )
    parser.add_argument(
        "-j", "--joint_prob_file", 
        help="A path to a temporally output file including a list of joint "\
             "probabilities of each mutation across sites."
    )
    parser.add_argument(
        "-p", "--position_file", 
        help="A path to position data file.", default=''
    )
    parser.add_argument(
        "-q", "--position_type", 
        help="Position type. Please specify 'table' or 'posmap'", 
        default=''
    )
    parser.add_argument(
        "-i", "--item_list_file", 
        help="A path to an item list file used for alignment concatenation in "\
             "BASEML.", 
        default=''
    )
    parser.add_argument(
        "-n", "--bootstrap_num", 
        help="Number of bootstrap resampling.", type=int, default=0
    )

    args = parser.parse_args()
    main(
        args.output_dir, 
        args.joint_prob_file, 
        args.position_file,
        args.position_type,
        args.item_list_file, 
        args.bootstrap_num
    )
