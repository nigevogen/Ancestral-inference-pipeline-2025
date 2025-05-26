open (IN1,"list_of_BTW_ancestor_prob_and_codon_path_prob.txt");
@str1 = <IN1>;

foreach $_ (@str1) {
    if ($_ =~ /(\d+)\t(\d+)\t(\d+)\t(\w\w\w)\t(\w\w\w)\t(.*)\t(.*)\t\t/) {
        $codon_pos = $1;
        $anc_node = $2;
        $der_node = $3;
        $anc_codon = $4;
        $der_codon = $5;
        $freq = $6;
        $prob_BTW = $7;
        
        if ($_ =~ /\((.*)\)\t(\d\.\d+)/) {
            $str = $1;
            $prob_CP = $2;
            $total_prob = $prob_BTW*$prob_CP;
            if ($str =~ /\'(\w\w\w)\'\,\s\'(\w\w\w)\'\,\s\'(\w\w\w)\'\,\s\'(\w\w\w)\'/) {
                open (OUT1, ">>list_of_codon_change_with_probability.txt");
                print (OUT1 "$codon_pos\t$anc_node\t$der_node\t$anc_codon\t$der_codon\t$freq\t$1\t$2\t$3\t$4\t$total_prob\n");
                close (OUT1);
            }
            elsif ($str =~ /\'(\w\w\w)\'\,\s\'(\w\w\w)\'\,\s\'(\w\w\w)\'/) {
                open (OUT1, ">>list_of_codon_change_with_probability.txt");
                print (OUT1 "$codon_pos\t$anc_node\t$der_node\t$anc_codon\t$der_codon\t$freq\t$1\t$2\t$3\t\t$total_prob\n");
                close (OUT1);
            }
            elsif ($str =~ /\'(\w\w\w)\'\,\s\'(\w\w\w)\'/) {
                open (OUT1, ">>list_of_codon_change_with_probability.txt");
                print (OUT1 "$codon_pos\t$anc_node\t$der_node\t$anc_codon\t$der_codon\t$freq\t$1\t$2\t\t\t$total_prob\n");
                close (OUT1);
            }
        }
    }
}
