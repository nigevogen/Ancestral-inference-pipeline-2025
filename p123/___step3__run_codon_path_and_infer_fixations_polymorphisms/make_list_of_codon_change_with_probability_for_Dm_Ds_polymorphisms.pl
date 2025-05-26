open (IN1,"list_of_codon_change_with_probability.txt");
@str1 = <IN1>;

foreach $_ (@str1) {
    if ($_ =~ /(\d+)\t7\t1\t(\w\w\w)\t(\w\w\w)\t(\d+)\t(\w\w\w.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $freq = $4;
        $line = $5;
        $prob = $6;
        
        #When a path contains 1 changes
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t\t$/) {
            $anc = $1;
            $der = $2;
            if ($anc ne $der) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_poly_r.txt");
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
                close (OUT1);
            }
        }
        #When a path contains >1 changes only the last change is considered as polymorphism
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t$/) {
            $anc = $2;
            $der = $3;
            if ($anc ne $der) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_poly_r.txt");
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
                close (OUT1);
                printf "a\n";
            }
        }
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)$/) {
            $anc = $4;
            $der = $5;
            if ($anc ne $der) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_poly_r.txt");
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
                close (OUT1);
                printf "b\n";
            }
        }
    }
    
    if ($_ =~ /(\d+)\t8\t2\t(\w\w\w)\t(\w\w\w)\t(\d+)\t(\w\w\w.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $freq = $4;
        $line = $5;
        $prob = $6;
        
        #When a path contains 1 changes
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t\t$/) {
            $anc = $1;
            $der = $2;
            if ($anc ne $der) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_poly_r.txt");
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
                close (OUT1);
            }
        }
        #When a path contains >1 changes only the last change is considered as polymorphism
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t$/) {
            $anc = $2;
            $der = $3;
            if ($anc ne $der) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_poly_r.txt");
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
                close (OUT1);
                printf "c\n";
            }
        }
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)$/) {
            $anc = $4;
            $der = $5;
            if ($anc ne $der) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_poly_r.txt");
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
                close (OUT1);
                printf "d\n";
            }
        }
    }
}
