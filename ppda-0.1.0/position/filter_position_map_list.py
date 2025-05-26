import argparse
from collections import namedtuple
from warnings import warn

def main(output_pos_map_file, input_pos_map_file, filt_pos_list_file, mode, 
         min_pos_num, apply_func=int):
    
    orig_pos_map = read_position_map_file(input_pos_map_file, apply_func)

    filt_pos_d = read_position_list_file(filt_pos_list_file, apply_func)

    filt_pos_map = filter_pos_map(orig_pos_map, filt_pos_d, mode, min_pos_num)

    save_as_position_map_file(output_pos_map_file, filt_pos_map)

def filter_pos_map(pos_map, filt_pos_d, mode, min_pos_num):
    if mode not in {'retain', 'remove'}:
        raise TypeError('Please specify "retain" or "remove" for filter mode. '\
                        f'{mode} is found.')
    new_keys = []
    filt_concat_pos_list = []
    filt_orig_pos_list = []

    for k, orig_pos_list in pos_map.orig_pos_d.items():
        tmp_concat_pos_list = []
        tmp_orig_pos_list = []

        try:
            tmp_filt_list = set(filt_pos_d[k])
        except KeyError:
            warn(f'item "{k}" is not found in the given position filter list.')
            continue

        for i, orig_pos in enumerate(orig_pos_list):
            if mode == 'retain':
                if orig_pos in tmp_filt_list:
                    tmp_concat_pos_list.append(pos_map.concat_pos_d[k][i])
                    tmp_orig_pos_list.append(orig_pos)
                    continue

            elif mode == 'remove':
                if orig_pos not in tmp_filt_list:
                    tmp_concat_pos_list.append(pos_map.concat_pos_d[k][i])
                    tmp_orig_pos_list.append(orig_pos)
                    continue
        
        if len(tmp_concat_pos_list) < min_pos_num:
            continue

        new_keys.append(k)
        filt_concat_pos_list.append(tmp_concat_pos_list)
        filt_orig_pos_list.append(tmp_orig_pos_list)

    return PositionMap(
        pos_map.prefix, 
        pos_map.suffix, 
        dict(zip(new_keys, filt_concat_pos_list)),
        dict(zip(new_keys, filt_orig_pos_list))
    )


def read_position_list_file(pos_list_file, apply_func=None):
    pos_d = {}
    tmp_pos_list = []
    item = ''

    with open(pos_list_file, 'r') as f:
        for l in f:
            if l[:-1] == '':
                continue

            if l.startswith('itemnum'):
                exp_itemnum = int(l[:-1].split('itemnum:')[1].strip())
                continue

            if l.startswith('>'):
                if item != '':
                    pos_d[apply_func(item)] = tmp_pos_list

                item = l[1:-1]
                tmp_pos_list = []
            else:
                tmp_pos_list.append(int(l[:-1]))

    if item != '':
        pos_d[apply_func(item)] = tmp_pos_list

    assert len(pos_d) == exp_itemnum, f'{len(pos_d)} != {exp_itemnum}'

    return pos_d

PositionMap = namedtuple(
    'PositionMap',
    ['prefix', 'suffix', 'concat_pos_d', 'orig_pos_d']
)

def read_position_map_file(pos_file, apply_func=None):
    if apply_func == None:
        apply_func = lambda s: s
        
    item_names = []
    concat_pos_list = []
    temp_concat_pos_list = []
    orig_pos_list = []
    temp_orig_pos_list = []
    
    with open(pos_file, 'r') as f:
        for l in f:
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
                if temp_concat_pos_list:
                    item_names.append(item)
                    concat_pos_list.append(temp_concat_pos_list)
                    orig_pos_list.append(temp_orig_pos_list)
                    
                item = apply_func(l[1:-1])
                temp_concat_pos_list = []
                temp_orig_pos_list = []
                
            else:
                temp_concat_pos_list.append(int(l[:-1].split('\t')[0]))
                temp_orig_pos_list.append(int(l[:-1].split('\t')[1]))
    
    if temp_concat_pos_list:
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

def save_as_position_map_file(output_file, pos_map):
    out_lines = []
    
    for item, concat_pos_list in pos_map.concat_pos_d.items():
        out_lines.append(
            '>{aln_name}\n{positions}'.format(
                aln_name=item,
                positions= '\n'.join([
                    f'{concat_pos}\t{pos_map.orig_pos_d[item][i]}'
                    for i, concat_pos in enumerate(concat_pos_list)
                ])
            )
        )
    with open(output_file, 'w') as f:
        print('prefix:', pos_map.prefix, file=f)
        print('suffix:', pos_map.suffix, file=f)
        print(f'itemnum: {len(pos_map.concat_pos_d)}', file=f)
        print('\n'.join(out_lines), file=f)

if __name__ == '__main__':
    desc = 'Filter position map file.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-o", "--output_pos_map_file", 
        help="A path to an output position map file."
    )
    parser.add_argument(
        "-i", "--input_pos_map_file", 
        help="A path to an input position map file."
    )
    parser.add_argument(
        "-f", "--filter_pos_file", 
        help="A path to a file of lists of positions to be retained/removed."
    )
    parser.add_argument(
        "-m", "--mode", 
        help="Mode of filtering, which controls whether given positions are "\
             "retained or removed in output position map. Please specify "\
             "'retain' or 'remove'."
    )
    parser.add_argument(
        "-n", "--min_pos_num", 
        help="Minimum number of sites from a single CDS to use divergence "\
             "calculation.", 
        type=int
    )

    args = parser.parse_args()
    main(
        args.output_pos_map_file, 
        args.input_pos_map_file,
        args.filter_pos_file, 
        args.mode, 
        args.min_pos_num,
        apply_func=int
    )
