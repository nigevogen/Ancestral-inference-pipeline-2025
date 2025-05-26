open (IN1,"joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW.txt");
@str1 = <IN1>;


$bookmark = 0;
if (@str1[$#str1] =~ /^(\d+)\t/) {
    $cod_num = $1;
}

for ($n1=0;$n1<=$cod_num;$n1++) {
    
    
    #read ancestral codon cofig file
    $prob_total = 0;
    $count_config = 0;
    
    for ($n2=$bookmark;$n2<=$#str1;$n2++) {
        
        if (@str1[$n2] =~ /^$n1\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\d.*)\n/) {
            
            $ms = $1;
            $ye = $2;
            $m = $3;
            $s = $4;
            $prob = $5;
            
            $check_ms = 0;
            $check_ye = 0;
            $check_m = 0;
            $check_s = 0;
            
            if (($ms eq 'TAA')||($ms eq 'TAG')||($ms eq 'TGA')) {
                $check_ms = 1;
            }
            if (($ye eq 'TAA')||($ye eq 'TAG')||($ye eq 'TGA')) {
                $check_ye = 1;
            }
            if (($m eq 'TAA')||($m eq 'TAG')||($m eq 'TGA')) {
                $check_m = 1;
            }
            if (($s eq 'TAA')||($s eq 'TAG')||($s eq 'TGA')) {
                $check_s = 1;
            }

            #When all paths are accessible by one or zero change
            if ($check_ms == 0) {
                if ($check_ye == 0) {
                    if ($check_m == 0) {
                        if ($check_s == 0) {
                            
                            ${'ms_'.$count_config} = $ms;
                            ${'ye_'.$count_config} = $ye;
                            ${'m_'.$count_config} = $m;
                            ${'s_'.$count_config} = $s;
                            ${'prob_'.$count_config} = $prob;
                            
                            $prob_total = $prob_total + $prob;
                            
                            $count_config++;
                        }
                    }
                }
            }
        
            #if the next line is the info of the next codon pos, output the info of the current codon pos
            if (@str1[$n2+1] =~ /^(\d+)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\d.*)\n/) {
                if ($1 == $n1+1) {
                    
                    if ($prob_total > 0) {
                        for ($n4=0;$n4<=$count_config-1;$n4++) {
                            $prob_control = ${'prob_'.$n4}/$prob_total;
                            
                            open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons.txt");
                            print (OUT1 "$n1\t${'ms_'.$n4}\t${'ye_'.$n4}\t${'m_'.$n4}\t${'s_'.$n4}\t$prob_control\n");
                            close (OUT1);
                        }
                    }
                    if ($prob_total == 0) {
                        open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons.txt");
                        print (OUT1 "$n1\tno_path\n");
                        close (OUT1);
                    }
                    $bookmark = $n2+1;
                    last;
                }
            }
            #if the next line is empty, output the info of the current codon pos
            if (@str1[$n2+1] !~ /^\d+/) {
                if ($prob_total > 0) {
                    for ($n4=0;$n4<=$count_config-1;$n4++) {
                        $prob_control = ${'prob_'.$n4};
                        
                        open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons.txt");
                        print (OUT1 "$n1\t${'ms_'.$n4}\t${'ye_'.$n4}\t${'m_'.$n4}\t${'s_'.$n4}\t$prob_control\n");
                        close (OUT1);
                    }
                }
                if ($prob_total == 0) {
                    open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons.txt");
                    print (OUT1 "$n1\tno_path\n");
                    close (OUT1);
                }
            }
        }
    }
}
