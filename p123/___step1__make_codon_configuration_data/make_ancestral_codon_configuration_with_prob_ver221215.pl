open (IN1,"data/anc_site_probs_all_node_1st_pos.txt");
@str1 = <IN1>;
open (IN2,"data/anc_site_probs_all_node_2nd_pos.txt");
@str2 = <IN2>;
open (IN3,"data/anc_site_probs_all_node_3rd_pos.txt");
@str3 = <IN3>;

for ($n=0;$n<=$#str1;$n++) {
    if (@str1[$n] =~ /^(\d+)\t/) {
        if (@str2[$n] =~ /^(\d+)\t/) {
            if (@str3[$n] =~ /^(\d+)\t/) {
                
                $cod_num = $1;
                
                $count_pos1 = 0;
                $count_pos2 = 0;
                $count_pos3 = 0;
                while (@str1[$n] =~ /(\w\t\w\t\w\t\w)\t(\d\.\d+)\t/g) {
                    ${'pos1_config'.$count_pos1} = $1;
                    ${'pos1_prob'.$count_pos1} = $2;
                    $count_pos1++;
                }
                while (@str1[$n] =~ /(\w\t\w\t\w\t\w)\t(\d\.\d+e\-\d+)\t/g) {
                    ${'pos1_config'.$count_pos1} = $1;
                    ${'pos1_prob'.$count_pos1} = $2;
                    $count_pos1++;
                }
                while (@str1[$n] =~ /(\w\t\w\t\w\t\w)\t(1)\t/g) {
                    ${'pos1_config'.$count_pos1} = $1;
                    ${'pos1_prob'.$count_pos1} = $2;
                    $count_pos1++;
                }
                
                while (@str2[$n] =~ /(\w\t\w\t\w\t\w)\t(\d\.\d+)\t/g) {
                    ${'pos2_config'.$count_pos2} = $1;
                    ${'pos2_prob'.$count_pos2} = $2;
                    $count_pos2++;
                }
                while (@str2[$n] =~ /(\w\t\w\t\w\t\w)\t(\d\.\d+e\-\d+)\t/g) {
                    ${'pos2_config'.$count_pos2} = $1;
                    ${'pos2_prob'.$count_pos2} = $2;
                    $count_pos2++;
                }
                while (@str2[$n] =~ /(\w\t\w\t\w\t\w)\t(1)\t/g) {
                    ${'pos2_config'.$count_pos2} = $1;
                    ${'pos2_prob'.$count_pos2} = $2;
                    $count_pos2++;
                }
                
                while (@str3[$n] =~ /(\w\t\w\t\w\t\w)\t(\d\.\d+)\t/g) {
                    ${'pos3_config'.$count_pos3} = $1;
                    ${'pos3_prob'.$count_pos3} = $2;
                    $count_pos3++;
                }
                while (@str3[$n] =~ /(\w\t\w\t\w\t\w)\t(\d\.\d+e\-\d+)\t/g) {
                    ${'pos3_config'.$count_pos3} = $1;
                    ${'pos3_prob'.$count_pos3} = $2;
                    $count_pos3++;
                }
                while (@str3[$n] =~ /(\w\t\w\t\w\t\w)\t(1)\t/g) {
                    ${'pos3_config'.$count_pos3} = $1;
                    ${'pos3_prob'.$count_pos3} = $2;
                    $count_pos3++;
                }
                #printf "$cod_num\t$count_pos1\t$count_pos2\t$count_pos3\n";
                
                for ($n1=0;$n1<=$count_pos1-1;$n1++) {
                    for ($n2=0;$n2<=$count_pos2-1;$n2++) {
                        for ($n3=0;$n3<=$count_pos3-1;$n3++) {
                            
                            if (${'pos1_config'.$n1} =~ /(\w)\t(\w)\t(\w)\t(\w)/) {
                                $pos1_node1 = $1;
                                $pos1_node2 = $2;
                                $pos1_node3 = $3;
                                $pos1_node4 = $4;
                            }
                            if (${'pos2_config'.$n2} =~ /(\w)\t(\w)\t(\w)\t(\w)/) {
                                $pos2_node1 = $1;
                                $pos2_node2 = $2;
                                $pos2_node3 = $3;
                                $pos2_node4 = $4;
                            }
                            if (${'pos3_config'.$n3} =~ /(\w)\t(\w)\t(\w)\t(\w)/) {
                                $pos3_node1 = $1;
                                $pos3_node2 = $2;
                                $pos3_node3 = $3;
                                $pos3_node4 = $4;
                            }
                            
                            $codon_node1 = "$pos1_node1$pos2_node1$pos3_node1";
                            $codon_node2 = "$pos1_node2$pos2_node2$pos3_node2";
                            $codon_node3 = "$pos1_node3$pos2_node3$pos3_node3";
                            $codon_node4 = "$pos1_node4$pos2_node4$pos3_node4";
                            $prob = ${'pos1_prob'.$n1}*${'pos2_prob'.$n2}*${'pos3_prob'.$n3};
                            
                            if ($prob > 0) {
                                open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW.txt");
                                print (OUT1 "$cod_num\t$codon_node1\t$codon_node2\t$codon_node3\t$codon_node4\t$prob\n");
                            }
                        }
                    }
                }
            }
        }
    }
}
