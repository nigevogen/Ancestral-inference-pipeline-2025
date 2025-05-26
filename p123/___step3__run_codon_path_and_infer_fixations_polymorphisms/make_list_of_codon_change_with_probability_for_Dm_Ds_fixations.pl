open (IN1,"list_of_codon_change_with_probability.txt");
@str1 = <IN1>;


foreach $_ (@str1) {
    #count changes occured in Dmel lineage
    if ($_ =~ /(\d+)\t5\t7\t(\w\w\w)\t(\w\w\w)\tbetween\sancestors\t(\w\w\w.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $line = $4;

        $prob = $5;

        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
        }
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
            if ($codon2 ne $codon3) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
                close (OUT1);
            }
        }
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            $codon4 = $4;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
            if ($codon2 ne $codon3) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
                close (OUT1);
            }
            if ($codon3 ne $codon4) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$codon3\t$codon4\t$prob\n");
                close (OUT1);
            }
        }
    }

    #count non-last changes occured between m' and mc1, mc2
    if ($_ =~ /(\d+)\t7\t1\t(\w\w\w)\t(\w\w\w)\t\d+\t(\w\w\w.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $line = $4;
        $prob = $5;

        #the last change is polymorphism and does not consider here
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
        }
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
            if ($codon2 ne $codon3) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage.txt");
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
                close (OUT1);
            }
        }
    }
    
    #count changes occured in Dsim lineage
    if ($_ =~ /(\d+)\t5\t8\t(\w\w\w)\t(\w\w\w)\tbetween\sancestors\t(\w\w\w.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $line = $4;
        $prob = $5;

        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
        }
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
            if ($codon2 ne $codon3) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
                close (OUT1);
            }
        }
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            $codon4 = $4;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
            if ($codon2 ne $codon3) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
                close (OUT1);
            }
            if ($codon3 ne $codon4) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$codon3\t$codon4\t$prob\n");
                close (OUT1);
            }
        }
    }

    #count non-last changes occured between m' and mc1, mc2
    if ($_ =~ /(\d+)\t8\t2\t(\w\w\w)\t(\w\w\w)\t\d+\t(\w\w\w.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $line = $4;
        $prob = $5;

        #the last change is polymorphism and does not consider here
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
        }
        if ($line =~ /^(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            if ($codon1 ne $codon2) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$$codon1\t$codon2\t$prob\n");
                close (OUT1);
            }
            if ($codon2 ne $codon3) {
                open (OUT1, ">>list_of_codon_change_with_probability_Dsim_lineage.txt");
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
                close (OUT1);
            }
        }
    }
}
