import os
import re
import sys
import numpy as np
from codonpaths import generate_codon_path_probs
from codonpaths import CODON_CHANGE_PATHS

# -------- INPUT ARGUMENTS -------- #
#specify species number
sp_num = 4
#specify ancestoral node number
anc_num = 4

# codonpaths parameter
synonymous_rate = 1.0
nonsynonymous_rate = 0.126

#give a list of derived..ancestor node pairs
#current setting is (1:Dmel, 2:Dsim, 3:Dyak, 4:Dere, 5:root, 6:ye, 7:m, 8:s)
ancestor_derived_pairs = ["1..7", "2..8", "3..6", "4..6", "6..5", "7..5", "8..5"]

# Paths
extant_path = 'codon_configuration_at_extant_nodes.txt'
ancestor_path = 'joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons.txt'
output_path = 'test_out.txt'

#get line number in ancestor file
# length = len(open('joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW.txt').readlines())

# ---------------------------------- #

if os.path.isfile(output_path):
    raise FileExistsError(f'Output file already exists: {output_path}.')

#from ancestor_dervied_pairs, make pairs of ancestor and derived codons 
# (anc_codon and der_codon), which will be used as the input of codon_path
target_node_list = [
    m for target_string in ancestor_derived_pairs 
        for m in target_string.split('..')
]

#run codon_path by KK using anc_codon and der_codon as the input
path_probs = generate_codon_path_probs(
    CODON_CHANGE_PATHS, synonymous_rate, nonsynonymous_rate,
    no_stop_codon_path=True)

#put the informations in ancestor file to the ancesotr_data
prev_codon_pos = -1
tmp_list = []
ancestor_data = {}

with open(ancestor_path, 'r') as ancestor:
    for line in ancestor:
        codon_pos = line.split('\t')[0]
        
        if prev_codon_pos != codon_pos:
            if tmp_list:
                ancestor_data[prev_codon_pos] = tmp_list
            tmp_list = []

        tmp_list.append(line)
        prev_codon_pos = codon_pos

ancestor_data[prev_codon_pos] = tmp_list

output_lines = []

#read extant file from the top
with open(extant_path, "r") as extant:
    for i, line in enumerate(extant):
        extant_items = line[:-1].split('\t')[:-1]

        assert len(extant_items[1:]) % 4 == 0, \
            f'Wrong formatting at extant line num {i}: {line}.'

        codon_pos = extant_items[0]

        if int(codon_pos) % 10000 == 0:
            print(codon_pos, end='\n')

        elif int(codon_pos) % 1000 == 0:
            print('.', end='')

        #put the informations in extant file to the extant_data
        extant_data_codon_1, extant_data_freq_1, extant_data_codon_2, extant_data_freq_2 = \
            zip(*[extant_items[i:i+4] for i in range(1, len(extant_items), 4)])
        
        #get the ancestor_data which has the same codon_pos index
        for anc_line in ancestor_data[codon_pos]:
            # if re.compile("^"+re.escape(codon_pos)+"\t").search(anc_line):

            ancestor_line = anc_line.split('\t')
            
            #put the informations in ancestor and extant codons in combined_data (make the codon configuration of extant+ancestors)
            #current setting is (m1, s1, y, e, root, ye, s, m) in combined_data_codon(freq)_1
            #               and (m2, s2, y, e, root, ye, s, m) in combined_data_codon(freq)_2
            combined_data_codon_1 = [
                extant_data_codon_1[n] if n < sp_num
                else ancestor_line[n-(sp_num-1)]
                    for n in range(sp_num+anc_num)
            ]
            combined_data_codon_2 = [
                extant_data_codon_2[n] if n < sp_num
                else ancestor_line[n-(sp_num-1)]
                    for n in range(sp_num+anc_num)
            ]
            combined_data_freq_1 = extant_data_freq_1
            combined_data_freq_2 = extant_data_freq_2

            #get the probability of the current codon configuration
            prob_n = ancestor_line[anc_num+1].rstrip("\n")
        
            n2 = 0
            tmp_output_lines = []
            while n2 < len(target_node_list):
                
                anc_node = int(target_node_list[n2+1])
                der_node = int(target_node_list[n2])
                anc_codon = combined_data_codon_1[anc_node-1]
                der_codon_1 = combined_data_codon_1[der_node-1]
                
                #when the derived node is extant, we consider the possibility of 
                # polymorphism (two codons)
                if (der_node<=sp_num):
                    der_codon_2 = combined_data_codon_2[der_node-1]
                    der_freq_1 = combined_data_freq_1[der_node-1]
                    der_freq_2 = combined_data_freq_2[der_node-1]
                n2 = n2 + 2
                
                # print(codon_pos, anc_node, anc_codon, 
                #     der_node, der_codon_1, prob_n+'\n')
                    
                if (der_node > sp_num):
                    codon_path_result = path_probs[(anc_codon, der_codon_1)]
                    codon_path_result_key_list = list(codon_path_result.keys())
                    codon_path_result_value_list = list(codon_path_result.values())

                    for n3 in range(len(codon_path_result_key_list)):
                        output_line = '\t'.join([str(a) for a in [
                            codon_pos, 
                            anc_node, der_node, 
                            anc_codon, der_codon_1, "between ancestors", prob_n, 
                            "", 
                            codon_path_result_key_list[n3], 
                            codon_path_result_value_list[n3]
                        ]])
                        tmp_output_lines.append(output_line)
                                                                    
                                                                    
                if (der_node <= sp_num):
                    codon_path_result = path_probs[(anc_codon, der_codon_1)]
                    codon_path_result_key_list = list(codon_path_result.keys())
                    codon_path_result_value_list = list(codon_path_result.values())

                    for n3 in range(len(codon_path_result_key_list)):
                        output_line = '\t'.join([str(a) for a in [
                            codon_pos, 
                            anc_node, der_node, 
                            anc_codon, der_codon_1, der_freq_1, prob_n, 
                            "", 
                            codon_path_result_key_list[n3], 
                            codon_path_result_value_list[n3]
                        ]])
                        tmp_output_lines.append(output_line)
                                                                                
                    if (der_codon_2 != "none"):
                        codon_path_result = path_probs[(anc_codon, der_codon_2)]
                        codon_path_result_key_list = list(codon_path_result.keys())
                        codon_path_result_value_list = list(codon_path_result.values())
                        
                        for n3 in range (0, len(codon_path_result_key_list)):
                            output_line = '\t'.join([str(a) for a in [
                                codon_pos, 
                                anc_node, der_node, 
                                anc_codon, der_codon_2, der_freq_2, prob_n, 
                                "", 
                                codon_path_result_key_list[n3], 
                                codon_path_result_value_list[n3]
                            ]])
                            tmp_output_lines.append(output_line)

            if int(codon_pos) % 10000 == 0:
                print('\n'.join(tmp_output_lines))

            output_lines += tmp_output_lines

with open(output_path, 'w') as f:
    print('\n'.join(output_lines), file=f)
