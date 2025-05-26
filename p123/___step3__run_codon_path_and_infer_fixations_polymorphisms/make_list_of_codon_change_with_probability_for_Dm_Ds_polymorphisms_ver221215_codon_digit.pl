open (IN1,"list_of_codon_change_with_probability_codon_digit.txt");
@str1 = <IN1>;

open (OUT1, ">>list_of_codon_change_with_probability_Dmel_poly_codon_digit.txt");
open (OUT2, ">>list_of_codon_change_with_probability_Dsim_poly_codon_digit.txt");

foreach $_ (@str1) {
    if ($_ =~ /(\d+)\t7\t1\t(\d+)\t(\d+)\t(\d+)\t(\d+.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $freq = $4;
        $line = $5;
        $prob = $6;
        
        #When a path contains 1 changes
        if ($line =~ /^(\d+)\t(\d+)\t\t$/) {
            $anc = $1;
            $der = $2;
            if ($anc ne $der) {
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
            }
        }
        #When a path contains >1 changes only the last change is considered as polymorphism
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t$/) {
            $anc = $2;
            $der = $3;
            if ($anc ne $der) {
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
            }
        }
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)$/) {
            $anc = $4;
            $der = $5;
            if ($anc ne $der) {
                print (OUT1 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
            }
        }
    }
    
    if ($_ =~ /(\d+)\t8\t2\t(\d+)\t(\d+)\t(\d+)\t(\d+.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $freq = $4;
        $line = $5;
        $prob = $6;
        
        #When a path contains 1 changes
        if ($line =~ /^(\d+)\t(\d+)\t\t$/) {
            $anc = $1;
            $der = $2;
            if ($anc ne $der) {
                print (OUT2 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
            }
        }
        #When a path contains >1 changes only the last change is considered as polymorphism
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t$/) {
            $anc = $2;
            $der = $3;
            if ($anc ne $der) {
                print (OUT2 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
            }
        }
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)$/) {
            $anc = $4;
            $der = $5;
            if ($anc ne $der) {
                print (OUT2 "$codon_pos\t$anc\t$der\t$freq\t$prob\n");
            }
        }
    }
}
