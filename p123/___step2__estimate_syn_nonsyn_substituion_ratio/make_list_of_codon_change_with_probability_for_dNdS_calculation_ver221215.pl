open (IN1,"joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_stop_codons_only_paths_with_single_or_no_change.txt");
@str1 = <IN1>;
open (IN2,"codon_configuration_at_extant_nodes.txt");
@str2 = <IN2>;

$bookmark = 0;
if (@str1[$#str1] =~ /^(\d+)\t/) {
    $cod_num = $1;
}

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
    
    for ($n2=$bookmark;$n2<=$#str1;$n2++) {
        
        if (@str1[$n2] =~ /^$n1\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\d.*)\n/) {
            $ms = $1;
            $ye = $2;
            $m = $3;
            $s = $4;
            $prob = $5;
            
            if ($ms ne $m) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage_for_NG_dNdS_calculation.txt");
                print (OUT1 "$n1\t$ms\t$m\t$prob\n");
                close (OUT1);
            }
            if ($ms ne $s) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage_for_NG_dNdS_calculation.txt");
                print (OUT1 "$n1\t$ms\t$s\t$prob\n");
                close (OUT1);
            }
            if ($ms ne $ye) {
                open (OUT1, ">>list_of_codon_change_with_probability_ye-ms_lineage_for_NG_dNdS_calculation.txt");
                print (OUT1 "$n1\t$ye\t$ms\t$prob\n");
                close (OUT1);
            }
            if ($ye ne $Dyak) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dyak_lineage_for_NG_dNdS_calculation.txt");
                print (OUT1 "$n1\t$ye\t$Dyak\t$prob\n");
                close (OUT1);
            }
            if ($ye ne $Dere) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dere_lineage_for_NG_dNdS_calculation.txt");
                print (OUT1 "$n1\t$ye\t$Dere\t$prob\n");
                close (OUT1);
            }
        }
        if (@str1[$n2+1] =~ /^(\d+)\t/) {
            if ($1 == $n1+1) {
                $bookmark = $n2+1;
                last;
            }
        }
    }
}

