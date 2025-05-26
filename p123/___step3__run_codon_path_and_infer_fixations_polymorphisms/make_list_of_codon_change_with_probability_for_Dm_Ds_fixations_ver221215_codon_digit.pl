open (IN1,"list_of_codon_change_with_probability_codon_digit.txt");
@str1 = <IN1>;

open (OUT1, ">>list_of_codon_change_with_probability_Dmel_lineage_codon_digit.txt");
open (OUT2, ">>list_of_codon_change_with_probability_Dsim_lineage_codon_digit.txt");

foreach $_ (@str1) {
    #count changes occured in Dmel lineage
    if ($_ =~ /(\d+)\t5\t7\t(\d+)\t(\d+)\tbetween\sancestors\t(\d+.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $line = $4;

        $prob = $5;

        if ($line =~ /^(\d+)\t(\d+)\t\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            if ($codon1 ne $codon2) {
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
        }
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            if ($codon1 ne $codon2) {
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
            if ($codon2 ne $codon3) {
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
            }
        }
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            $codon4 = $4;
            if ($codon1 ne $codon2) {
                print (OUT1 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
            if ($codon2 ne $codon3) {
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
            }
            if ($codon3 ne $codon4) {
                print (OUT1 "$codon_pos\t$codon3\t$codon4\t$prob\n");
            }
        }
    }

    #count non-last changes occured between m' and mc1, mc2
    if ($_ =~ /(\d+)\t7\t1\t(\d+)\t(\d+)\t\d+\t(\d+.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $line = $4;
        $prob = $5;

        #the last change is polymorphism and does not consider here
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            printf "$codon_pos a\n";
            if ($codon1 ne $codon2) {
                print (OUT1 "$codon_pos\t$$codon1\t$codon2\t$prob\n");
            }
        }
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            printf "$codon_pos b\n";
            if ($codon1 ne $codon2) {
                print (OUT1 "$codon_pos\t$$codon1\t$codon2\t$prob\n");
            }
            if ($codon2 ne $codon3) {
                print (OUT1 "$codon_pos\t$codon2\t$codon3\t$prob\n");
            }
        }
    }
    
    #count changes occured in Dsim lineage
    if ($_ =~ /(\d+)\t5\t8\t(\d+)\t(\d+)\tbetween\sancestors\t(\d+.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $line = $4;
        $prob = $5;

        if ($line =~ /^(\d+)\t(\d+)\t\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            if ($codon1 ne $codon2) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
        }
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            if ($codon1 ne $codon2) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
            if ($codon2 ne $codon3) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
        }
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            $codon4 = $4;
            if ($codon1 ne $codon2) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
            if ($codon2 ne $codon3) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
            if ($codon3 ne $codon4) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
        }
    }

    #count non-last changes occured between s' and sc1, sc2
    if ($_ =~ /(\d+)\t8\t2\t(\d+)\t(\d+)\t\d+\t(\d+.*)\t(\d.*)\n/) {
        $codon_pos = $1;
        $line = $4;
        $prob = $5;

        #the last change is polymorphism and does not consider here
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t$/) {
            $codon1 = $1;
            $codon2 = $2;
            printf "$codon_pos c\n";
            if ($codon1 ne $codon2) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
        }
        if ($line =~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)$/) {
            $codon1 = $1;
            $codon2 = $2;
            $codon3 = $3;
            printf "$codon_pos d\n";
            if ($codon1 ne $codon2) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
            if ($codon2 ne $codon3) {
                print (OUT2 "$codon_pos\t$codon1\t$codon2\t$prob\n");
            }
        }
    }
}
