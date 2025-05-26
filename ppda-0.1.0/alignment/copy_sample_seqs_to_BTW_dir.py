
import os
import argparse
from collections import namedtuple
from shutil import copy2

def main(orig_dir, btw_dir):
    # Copy alignments
    item_list_file = os.path.join(orig_dir, '0.itemlist')
    itemlist = read_item_list(item_list_file)
    file_names = reconstruct_file_names(itemlist)
    seq_dir = os.path.join(btw_dir, '_seqs_folder')
    out_flist = copy_alignments(orig_dir, file_names, seq_dir)
    print(len(out_flist), f'files were saved in {seq_dir}.')

    # Make binlist for BASEML
    bin_list_file = os.path.join(btw_dir, 'bin_dat_1_index_lists')
    to_item_list(
        bin_list_file, 'concate_seq', itemlist.item_list, 
        itemlist.prefix, itemlist.suffix
    )

    # Make list for concatenation in BTW
    list_path = os.path.join(btw_dir, '_seq_list_to_be_concatenated/seq_list')
    output_item_list_for_concatenation(list_path, file_names)

def copy_alignments(orig_dir, file_names, seq_dir):
    if not os.path.isdir(seq_dir):
        raise FileNotFoundError(seq_dir)

    for fname in file_names:
        orig_path = os.path.join(orig_dir, fname)
        dist_path = os.path.join(seq_dir, fname)

        if os.path.isfile(dist_path):
            raise FileExistsError(dist_path)

        copy2(orig_path, dist_path)

    flist = to_filelist(seq_dir)
    return flist

def reconstruct_file_names(itemlist):
    return [
        f'{itemlist.prefix}{item}{itemlist.suffix}'
        for item in itemlist.item_list
    ]

def output_item_list_for_concatenation(output_path, file_list):
    if os.path.isfile(output_path):
        raise FileExistsError(output_path)

    with open(output_path, 'w') as f:
        print('\n'.join(file_list), file=f)

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
    if os.path.isfile(out_path):
        raise FileExistsError(out_path)

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

def to_filelist(dir_path):
    """ Returns a file list in a given directory.  """
    l1 = os.listdir(dir_path)
    l2= [a for a in l1 if not a.startswith('.')]
    flist = [a for a in l2 if not a.startswith('0')]

    with open(os.path.join(dir_path, '0.filelist'), 'w') as f:
        print('itemnum: '+str(len(flist)), file=f)
        print('\n'.join(flist), file=f)
    return flist

if __name__ == '__main__':
    desc = 'Copy alignments and make concatenation file list in BTW directory.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "aln_dir", 
        help='A path to an alignment directory. Sequence alignments listed in '\
             '"0.itemlist" in this folder will be copied to <BTW_DIR>/_seqs_folder.')
    parser.add_argument(
        "btw_dir", 
        help='A path to BTW directory, which must include directories of '\
             '"_seqs_folder" and "_seq_list_to_be_concatenated".')

    # Parse arguments from terminal
    args = parser.parse_args()
    
    main(args.aln_dir, args.btw_dir)
