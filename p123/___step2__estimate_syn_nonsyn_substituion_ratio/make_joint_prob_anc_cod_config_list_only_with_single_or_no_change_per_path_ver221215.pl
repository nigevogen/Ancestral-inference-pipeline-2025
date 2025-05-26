open (IN1,"joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons.txt");
@str1 = <IN1>;
open (IN2,"codon_configuration_at_extant_nodes.txt");
@str2 = <IN2>;

$bookmark = 0;
if (@str1[$#str1] =~ /^(\d+)\t/) {
    $cod_num = $1;
}
printf "$cod_num\n";
for ($n1=0;$n1<=$cod_num;$n1++) {
    
    #get information of Dmel Dsim polymorphic states
    $extant_nodes = @str2[$n1];
    @extant_nodes = split /\t/, $extant_nodes;
    $Dmel_1 = @extant_nodes[1];
    $Dmel_2 = @extant_nodes[3];
    $Dsim_1 = @extant_nodes[5];
    $Dsim_2 = @extant_nodes[7];
    $Dyak = @extant_nodes[9];
    $Dere = @extant_nodes[13];
    
    @Dmel_1 = split //, $Dmel_1;
    @Dmel_2 = split //, $Dmel_2;
    @Dsim_1 = split //, $Dsim_1;
    @Dsim_2 = split //, $Dsim_2;
    @Dyak = split //, $Dyak;
    @Dere = split //, $Dere;
    
    #printf "@Dmel_1\t@Dmel_2\t@Dsim_1\t@Dsim_2\t@Dyak\t@Dere\n";
    
    
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
            
            @ms = split //, $ms;
            @ye = split //, $ye;
            @m = split //, $m;
            @s = split //, $s;
            
            $check_msye = 0;
            $check_msm = 0;
            $check_mss = 0;
            $check_yeDyak = 0;
            $check_yeDere = 0;
            $check_mDmel_1 = 0;
            $check_mDmel_2 = 0;
            $check_sDsim_1 = 0;
            $check_sDsim_2 = 0;
            
            #count #changes for each adjascent nodes
            for ($n3=0;$n3<=2;$n3++) {
                if (@ms[$n3] ne @ye[$n3]) {
                    $check_msye++;
                }
                if (@ms[$n3] ne @m[$n3]) {
                    $check_msm++;
                }
                if (@ms[$n3] ne @s[$n3]) {
                    $check_mss++;
                }
                if (@ye[$n3] ne @Dyak[$n3]) {
                    $check_yeDyak++;
                }
                if (@ye[$n3] ne @Dere[$n3]) {
                    $check_yeDere++;
                }
                if (@m[$n3] ne @Dmel_1[$n3]) {
                    $check_mDmel_1++;
                }
                if ($Dmel_2 ne 'none') {
                    if (@m[$n3] ne @Dmel_2[$n3]) {
                        $check_mDmel_2++;
                    }
                }
                if (@s[$n3] ne @Dsim_1[$n3]) {
                    $check_sDsim_1++;
                }
                if ($Dsim_2 ne 'none') {
                    if (@s[$n3] ne @Dsim_2[$n3]) {
                        $check_sDsim_2++;
                    }
                }
            }
           # printf "$check_msye\t$check_msm\t$check_mss\t$check_yeDyak\t$check_yeDere\t$check_mDmel_1\t$check_mDmel_1\t$check_sDsim_1\t$check_sDsim_2\n";
            #When all paths are accessible by one or zero change
            if ($check_msye <= 1) {
                if ($check_msm <= 1) {
                    if ($check_mss <= 1) {
                        if ($check_yeDyak <= 1) {
                            if ($check_yeDere <= 1) {
                                if ($check_mDmel_1 <= 1) {
                                    if ($check_mDmel_2 <= 1) {
                                        if ($check_sDsim_1 <= 1) {
                                            if ($check_sDsim_2 <= 1) {
                                                
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
                            }
                        }
                    }
                }
            }
            #if the next line is the info of the next codon pos, output the info of the current codon pos
            if (@str1[$n2+1] =~ /^(\d+)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\d.*)\n/) {
                if ($1 == $n1+1) {
                    
                    if ($prob_total > 0) {
                        for ($n4=0;$n4<=$count_config-1;$n4++) {
                            $prob_control = ${'prob_'.$n4};
                            
                            open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_stop_codons_only_paths_with_single_or_no_change.txt");
                            print (OUT1 "$n1\t${'ms_'.$n4}\t${'ye_'.$n4}\t${'m_'.$n4}\t${'s_'.$n4}\t$prob_control\n");
                            close (OUT1);
                        }
                    }
                    if ($prob_total == 0) {
                        open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_stop_codons_only_paths_with_single_or_no_change.txt");
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
                        
                        open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_stop_codons_only_paths_with_single_or_no_change.txt");
                        print (OUT1 "$n1\t${'ms_'.$n4}\t${'ye_'.$n4}\t${'m_'.$n4}\t${'s_'.$n4}\t$prob_control\n");
                        close (OUT1);
                    }
                }
                if ($prob_total == 0) {
                    open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_stop_codons_only_paths_with_single_or_no_change.txt");
                    print (OUT1 "$n1\tno_path\n");
                    close (OUT1);
                }
            }
        }
    }
}
