import re
import argparse
import pandas as pd
from datetime import datetime
from collections import namedtuple, defaultdict

def main(output_file, anc_prob_file, mlf_file, joint_prob_file, internal_node_num, 
         poly_info_str, pos_list_file=''):
    """ Output summary of ancestor states. 

    Example
    -------
    /* Created by output_AI_result_summary.py on 2021-01-08T11:44:24.872929.
    /* 
    /* The following arguments are given,
    /* 	output_file = /Users/haruka/Desktop/test.txt
    /* 	anc_prob_file = /Volumes/1TB_ssd_SDa/Dropbox/Documents_DB/01_Projects/1_Data_Analysis/2_PopGen_analysis/5_Post_AI_site_choice_analysis_mpspye_r624/1_Analysis/210106_pre_iAA/1_2_0_2_Filtered_X_base_iAA/G/3_btw_result/anc_site_probs_all_node_2nd_BTW_1.txt
    /* 	mlf_file = /Volumes/1TB_ssd_SDa/Dropbox/Documents_DB/01_Projects/1_Data_Analysis/2_PopGen_analysis/5_Post_AI_site_choice_analysis_mpspye_r624/1_Analysis/210106_pre_iAA/1_2_0_2_Filtered_X_base_iAA/G/2_baseml_result/04_output/sample_result/result_bin_1/mlf_bin1
    /* 	joint_prob_file = /Volumes/1TB_ssd_SDa/Dropbox/Documents_DB/01_Projects/1_Data_Analysis/2_PopGen_analysis/5_Post_AI_site_choice_analysis_mpspye_r624/1_Analysis/210106_pre_iAA/1_2_0_2_Filtered_X_base_iAA/G/4_AWP_mutation_count/join_prob_with_freq.txt
    /* 	internal_node_num = 4
    /* 	poly_info_str = RG:RG_collapse_:3,4:9:0:14..MD:MD_collapse_:5,6:10:14:21
    /* 
    total_seq_len: 40061
    total_var_site_num: 24886
    total_extant_config_num: 979
    tree_from_mlf: [[7, 8], [8, 1], [8, 2], [7, 9], [9, 3], [9, 4], [7, 10], [10, 5], [10, 6]]

    base_composition:
        node 1:	T=18.42	C=52.49	A=22.52	G=6.56
        node 2:	T=17.1	C=51.81	A=23.65	G=7.45
        node 3:	T=19.92	C=49.22	A=24.38	G=6.48
        node 4:	T=19.91	C=49.19	A=24.41	G=6.49
        node 5:	T=18.6	C=49.87	A=24.18	G=7.36
        node 6:	T=18.53	C=49.98	A=24.16	G=7.33
        node 7:	T=17.4	C=51.63	A=23.95	G=7.02
        node 8:	T=17.05	C=51.27	A=24.3	G=7.38
        node 9:	T=19.09	C=50.18	A=24.29	G=6.44
        node 10:	T=17.75	C=51.28	A=24.0	G=6.97

    polymorphism_info:
        PopulationSampleInfo(sample_prefix='RG', collapse_prefix='RG_collapse_', extant_nodes=[3, 4], ancestor_node=9, sample_seq_start=0, allele_count=14)
        PopulationSampleInfo(sample_prefix='MD', collapse_prefix='MD_collapse_', extant_nodes=[5, 6], ancestor_node=10, sample_seq_start=14, allele_count=21)

    total_fixed_mutations:
        RG:	TC=90.221,	TA=56.468,	TG=16.697,	CT=525.864,	CA=209.334,	CG=75.709,	AT=170.746,	AC=83.595,	AG=57.537,	GT=142.921,	GC=54.215,	GA=183.921
        MD:	TC=153.483,	TA=49.884,	TG=19.127,	CT=243.567,	CA=121.87,	CG=50.912,	AT=77.691,	AC=78.611,	AG=73.916,	GT=40.026,	GC=43.618,	GA=77.994

    total_polymorphic_mutations:
        RG:	TC=77.346,	TA=28.085,	TG=12.236,	CT=349.654,	CA=106.318,	CG=49.051,	AT=73.915,	AC=25.682,	AG=56.995,	GT=26.264,	GC=13.949,	GA=60.975
        MD:	TC=276.275,	TA=93.61,	TG=57.024,	CT=575.725,	CA=260.252,	CG=135.32,	AT=134.39,	AC=103.248,	AG=176.712,	GT=44.476,	GC=48.18,	GA=128.288

    frequency_distribution:
        RG:
            TC:	41.15	9.428	3.119	1.363	3.099	3.21	1.581	0.866	1.495	3.551	1.721	1.263	5.501
            TA:	12.129	6.141	1.98	1.085	0.026	nan	0.742	0.899	0.439	1.88	1.208	nan	1.555
            TG:	6.019	1.499	nan	0.413	0.498	nan	nan	2.17	0.5	nan	0.695	nan	0.441
            CT:	172.999	45.237	20.779	23.449	15.005	9.634	9.919	7.79	6.401	9.137	6.381	8.072	14.85
            CA:	60.375	11.234	7.495	5.624	4.392	5.496	1.988	3.083	0.927	0.999	2.309	0.464	1.932
            CG:	32.207	5.412	2.0	1.5	1.472	0.318	1.0	0.977	nan	0.671	0.883	1.796	0.814
            AT:	27.945	12.5	9.292	3.62	1.561	3.101	1.758	1.0	1.474	2.915	2.02	2.859	3.871
            AC:	10.568	1.536	2.691	1.501	2.573	0.417	0.012	0.504	1.608	0.876	0.505	0.766	2.125
            AG:	22.735	5.532	2.448	0.964	1.076	1.227	2.582	1.626	1.407	0.405	0.94	nan	16.054
            GT:	9.559	5.0	1.805	1.5	0.5	0.83	nan	nan	1.002	1.587	3.0	1.001	0.481
            GC:	4.186	3.204	2.617	0.329	1.0	0.523	nan	0.682	1.028	nan	nan	0.088	0.293
            GA:	13.446	10.5	4.53	7.095	2.093	1.374	1.418	1.273	2.424	2.536	4.052	2.468	7.765
        MD:
            TC:	117.506	38.595	15.497	12.905	9.045	9.014	5.445	7.934	5.178	3.791	4.603	4.237	8.552	4.146	2.942	5.417	3.346	5.613	3.264	9.244
            TA:	49.041	15.264	6.076	2.767	3.5	1.899	0.499	0.82	0.437	1.515	0.886	nan	1.008	1.248	1.689	1.429	0.968	0.899	1.749	1.914
            TG:	30.499	8.032	3.995	3.0	1.476	0.916	nan	0.773	1.0	1.0	1.0	1.0	1.0	nan	nan	1.006	nan	0.97	0.209	1.149
            CT:	375.256	81.736	34.387	14.154	13.083	7.558	4.854	0.948	6.263	3.397	3.209	2.822	2.066	5.055	2.486	2.455	3.595	3.503	1.905	6.994
            CA:	172.47	32.402	15.939	5.143	7.681	4.404	1.648	3.0	2.974	1.489	2.113	1.913	0.276	1.157	2.208	1.202	1.754	0.354	0.702	1.422
            CG:	101.5	12.036	3.989	2.499	3.5	2.28	0.974	1.0	1.495	nan	0.494	0.601	nan	0.411	1.915	0.0	nan	0.762	0.995	0.871
            AT:	72.086	19.751	7.101	6.032	3.071	2.811	4.752	0.992	1.5	1.114	0.985	2.063	1.68	1.501	2.101	nan	2.233	1.424	2.736	0.459
            AC:	47.578	11.298	12.146	5.246	2.798	2.792	1.843	0.724	0.587	1.887	0.511	2.026	nan	0.852	1.096	2.319	2.357	2.061	2.598	2.53
            AG:	107.836	20.702	10.677	8.004	2.855	3.494	3.516	2.262	2.062	2.0	0.972	2.996	0.526	0.997	0.716	0.275	1.727	0.374	2.124	2.598
            GT:	24.351	7.791	4.03	nan	2.494	2.0	1.0	nan	nan	nan	0.5	nan	0.727	nan	0.584	0.524	nan	0.005	0.468	0.001
            GC:	23.129	6.505	3.738	2.5	2.5	1.085	2.089	1.0	0.399	0.506	nan	0.505	nan	2.026	0.72	nan	0.501	0.511	0.464	nan
            GA:	70.402	24.376	9.126	5.273	3.225	3.284	0.503	1.974	0.504	2.028	nan	1.938	1.238	0.484	1.506	0.645	0.496	0.323	0.798	0.164

    extant_configurations:
        AAAAAA: 5848
        AAAAAC: 63
        AAAAAG: 87
        ...
        
    """
    start = datetime.now().isoformat()
    args = {
        'output_file': output_file, 
        'anc_prob_file': anc_prob_file, 
        'mlf_file': mlf_file, 
        'joint_prob_file': joint_prob_file, 
        'internal_node_num': internal_node_num, 
        'poly_info_str': poly_info_str,
        'pos_list_file': pos_list_file
    }
    desc = make_description_line(start, **args)

    if pos_list_file != '':
        concat_pos = get_concat_pos_list(pos_list_file)
        concat_pos_set = set(concat_pos)
        assert len(concat_pos) == len(concat_pos_set)
    else:
        concat_pos_set = set()

    poly_info_list = poly_info_parser_from_terminal(poly_info_str)
    anc_info = read_anc_site_probs(
        anc_prob_file, len(poly_info_list), internal_node_num, concat_pos_set)
    base_compos = get_base_composition(anc_info)
    
    tree = get_branch_info_from_mlf(mlf_file)

    df = read_joint_prob_file(joint_prob_file, concat_pos_set)
    sfs_list = get_sfs_d_for_each_pop(df, poly_info_list)

    save_to_file(
        output_file, anc_info, base_compos, tree, poly_info_list, 
        sfs_list, desc
    )

def make_description_line(start_date, **kwargs):
    lines = [
        f'Created by output_AI_result_summary.py on {start_date}.\n',
        'The following arguments are given,',
    ]
    lines += [f'\t{k} = {str(v)}' for k, v in kwargs.items()]
    lines.append('')

    return '/* '+'\n/* '.join('\n'.join(lines).split('\n'))

def save_to_file(output_file, anc_info, base_compos, tree, poly_info_list, 
                 sfs_list, desc):
    output_lines = [
        desc,
        f'total_seq_len: {anc_info.total_seq_len}', 
        f'total_var_site_num: {anc_info.total_var_site_num}', 
        f'total_extant_config_num: {anc_info.total_extant_config_num}', 
        f'tree_from_mlf: {tree}', '', 'base_composition:'
    ]
    bases = 'TCAG'
    mutations = [b1+b2 for b1 in bases for b2 in bases if b1 != b2]
    mutations_t = sfs_list[0].index
    nodes = sorted(set([n for b in tree for n in b]))
    
    for node in nodes:
        dat = [
            '{}={}'.format(base, myround(base_compos[(node, base)], 2)) 
            for base in bases
        ]
        output_lines.append('\tnode {}:\t{}'.format(node, '\t'.join(dat)))
    
    fix_list = []
    poly_list = []
    fd_list = []

    output_lines.append(
        '\npolymorphism_info:\n\t{}'.format('\n\t'.join([
            str(poly_info)
            for poly_info in poly_info_list
        ]))
    )

    for poly_info, table in zip(poly_info_list, sfs_list):
        fix = f'\t{poly_info.sample_prefix}:\t' \
            + ',\t'.join([
            ''.join(t)+f'={myround(count, 3)}' 
            for t, count in table.loc[mutations_t, poly_info.allele_count].items()
        ])
        poly = f'\t{poly_info.sample_prefix}:\t' \
            + ',\t'.join([
            ''.join(t)+f'={myround(count, 3)}' 
            for t , count in table.loc[mutations_t, range(1, poly_info.allele_count)]\
                .sum(axis=1).items()
        ])
        fd = f'\t{poly_info.sample_prefix}:\n' \
            + '\n'.join([
            '\t\t{}:\t{}'.format(m, '\t'.join([str(myround(a, 3)) for a in dist])) 
            for m, dist in zip(
                mutations, 
                table.loc[mutations_t, range(1, poly_info.allele_count)].values)
        ])
        fix_list.append(fix)
        poly_list.append(poly)
        fd_list.append(fd)

    output_lines.append('\ntotal_fixed_mutations:')
    for fix in fix_list:
        output_lines.append(fix)

    output_lines.append('\ntotal_polymorphic_mutations:')
    for poly in poly_list:
        output_lines.append(poly)

    output_lines.append('\nfrequency_distribution:')
    for fd in fd_list:
        output_lines.append(fd)

    out = '\n'.join(sorted([
        f'\t{ex}: {len(pos)}' for ex, pos in anc_info.config_pos_d.items()]))
    output_lines.append(f'\nextant_configurations:\n{out}')

    with open(output_file, 'w') as f:
        print('\n'.join(output_lines), file=f)

def get_sfs_d_for_each_pop(df, poly_info_list):
    sfs_list = []
    pivot_table_kw = {
        'values': ['count'], 'index': ['anc_state', 'der_state'], 
        'columns': ['frequency']
    }

    for poly_info in poly_info_list:
        dat = df[
            (df['anc_node'] == poly_info.ancestor_node) |
            (df['der_node'] == poly_info.ancestor_node) 
        ].groupby(by=['anc_state', 'der_state', 'frequency']).sum()

        sfs_list.append(dat.pivot_table(**pivot_table_kw)['count'].fillna(0))

    return sfs_list

def get_base_composition(anc_info):
    """ Returns a dictionary of base composition for each node.
    Example
    -------
    {
        (1, T): 20, # (node_num, base): percentage across sites
        (1, C): 20, 
        (1, A): 30, 
        (1, G): 30, 
        (2, T): 23, 
        (2, C): 34, 
        ...
    }
    """
    base_compos_d = defaultdict(int)

    for ext, pos_list in anc_info.config_pos_d.items():
        ext_count = len(ext)
        site_num = len(pos_list)

        for i, base in enumerate(ext):
            base_compos_d[(i+1, base)] += site_num

        anc_probs = anc_info.config_anc_prob_d[ext]
        for internal, prob in anc_probs:
            for i, base in enumerate(internal):
                try:
                    base_compos_d[(ext_count+i+1, base)] += prob
                except:
                    raise TypeError(anc_probs)
                
    return {
        t: count / anc_info.total_seq_len * 100
        for t, count in base_compos_d.items()
    }

def get_branch_info_from_mlf(mlf_file):
    matches = []
    with open(mlf_file, 'r') as f:
        for l in f:
            if re.match(r'^\s+(\d+\.\.\d+\s+)+$', l):
               matches.append(l.rstrip()) 
    assert len(matches) == 1
    match = matches[0]

    return [[int(c) for c in a.split('..')] for a in match.split()]

AncestorStatesInfo = namedtuple(
    'AncestorStatesInfo',
    [
        'total_seq_len', # int
        'total_var_site_num', # int
        'total_extant_config_num', # int
        'config_pos_d', # dict
        'config_anc_prob_d', # dict. {'CCCCCT': [('CCCC', 0.999216467993116), ('CCCT', 0.000783532006883675)]}
    ]
)

def read_anc_site_probs(file_path, poly_sp_num, internal_node_num, concat_pos):
    """ Reads an output file of BTW. 
    """
    pos_d = defaultdict(list)
    prob_d = defaultdict(list)
    total_site_count = 0
    total_var_site_num = 0

    with open(file_path, 'r') as f:
        for l in f:
            pos, ex, anc = read_ancestor_state(l[:-1], poly_sp_num, internal_node_num)

            if len(concat_pos) > 0:
                if pos not in concat_pos:
                    continue

            pos_d[ex].append(pos)
            prob_d[ex] += anc

            total_site_count += 1

            if len(set(list(ex))) > 1:
                total_var_site_num += 1

    config_num = len(pos_d)
    assert len(prob_d) == config_num
    
    return AncestorStatesInfo(
        total_site_count, total_var_site_num, config_num, pos_d, prob_d)

def read_ancestor_state(s, poly_sp_num, internal_node_num):
    """ Reads a line of ancestor probability. """
    parts = s.split()
    pos = int(parts[0])

    poly_info_end_ix = int(poly_sp_num*2)+1
    extant = parts[poly_info_end_ix][:-1]
    anc = parts[poly_info_end_ix+1:]
    n = internal_node_num+1

    to_anc_prob = lambda s: (''.join(s[:internal_node_num]), float(s[internal_node_num]))

    anc_probs = [
        to_anc_prob(anc[n*i:n*i+n]) 
            for i in range(int(len(anc) / n))
    ]
    
    return pos, extant, anc_probs

def read_joint_prob_file(file_path, concat_pos):
    pos_list = []
    node1_list = []
    node2_list = []
    state1_list = []
    state2_list = []
    freq_list = []
    prob_list = []
    changes = []
    pos = None
    
    with open(file_path, 'r') as f:
        for l in f:
            if l.startswith('>'):
                if changes:
                    for n1, n2, m, freq, p in changes:
                        pos_list.append(pos)
                        node1_list.append(n1)
                        node2_list.append(n2)
                        state1_list.append(m[0])
                        state2_list.append(m[1])
                        freq_list.append(freq)
                        prob_list.append(p)
                        
                pos = int(l[1:-1])

            elif len(concat_pos) > 0:
                if pos in concat_pos:
                    changes = parse_probs(l[:-1])
                else:
                    changes = []

            else:
                changes = parse_probs(l[:-1])

    if changes:
        for n1, n2, m, freq, p in changes:
            pos_list.append(pos)
            node1_list.append(n1)
            node2_list.append(n2)
            state1_list.append(m[0])
            state2_list.append(m[1])
            freq_list.append(freq)
            prob_list.append(p)
    
    df = pd.DataFrame(
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
    cols = ['anc_node', 'der_node', 'anc_state', 'der_state', 'frequency']
    return df\
        .groupby(cols)['prob']\
        .sum().reset_index().rename(columns={'prob': 'count'})

def parse_probs(s):
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

PopulationSampleInfo = namedtuple(
    'PopulationSampleInfo', [
        'sample_prefix', # Prefix in sample seq
        'collapse_prefix', # Prefix in collapse seq
        'extant_nodes', # Order among collapse sequneces
        'ancestor_node', # Ancestor node id in a tree in mlf file (BASEML output)
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

def get_concat_pos_list(pos_file):
    item_names = []
    concat_pos_list = []
    item = ''
    
    with open(pos_file, 'r') as f:
        for l in f:
            if l[:-1] == '':
                continue

            if l.startswith('prefix'):
                continue
                
            if l.startswith('suffix'):
                continue
                
            if l.startswith('itemnum'):
                exp_itemnum = int(l[:-1].split('itemnum: ')[1])
                continue
                
            if l.startswith('>'):
                if item != '':
                    item_names.append(item)

                item = l[1:-1]
            else:
                pos = l[:-1].split('\t')
                assert len(pos) == 2, pos
                concat_pos_list.append(int(pos[0]))
    
    if item != '':
        item_names.append(item)

    assert len(item_names) == exp_itemnum, f'{len(item_names)} != {exp_itemnum}'
    
    return concat_pos_list

def myround(a, ndigits=2):
    n = 10 ** ndigits
    return (a * n * 2 + 1) // 2 / n

if __name__ == '__main__':
    desc = 'Output a file of Ancestor Inference (AI) summary. This script reads '\
        'outputs from BASEML (to read tree config), BTW (to read ancestor '\
        'states) and joint probability list (to calculate SFS).'

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-o", "--output_file", 
        help="A path to output file for AI result summary."
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
        "-j", "--joint_prob_file", 
        help="A path to an output file including a list of joint "\
             "probabilities of each mutation across sites."
    )
    parser.add_argument(
        "-i", "--internal_node_num", 
        help="Number of internal nodes in a sample tree used in BASEML run.",
        type=int
    )
    parser.add_argument(
        "-p", "--poly_info", 
        help="Information for population samples. Example: RG:RG_collapse_"\
             ":3,4:9:0:14..MD:MD_collapse_:5,6:10:14:21"
    )
    parser.add_argument(
        "-l", "--pos_list_file", 
        help="A path to position list file.",
        nargs='?', const='', default=''
    )
    args = parser.parse_args()
    main(
        args.output_file,
        args.anc_prob_file, 
        args.mlf_file,
        args.joint_prob_file, 
        args.internal_node_num, 
        args.poly_info, 
        args.pos_list_file
    )