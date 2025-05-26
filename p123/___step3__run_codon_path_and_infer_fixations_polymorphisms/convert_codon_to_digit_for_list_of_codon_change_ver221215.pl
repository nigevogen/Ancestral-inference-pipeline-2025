open (IN1,"list_of_codon_change_with_probability.txt");
@str1 = <IN1>;

open (OUT1, ">>list_of_codon_change_with_probability_codon_digit.txt");

foreach $_ (@str1) {
    @line = split /\t/, $_;
    
    for ($n=0;$n<=$#line;$n++) {
        if (@line[$n] eq 'TTT') {@line[$n] = 1;}
        if (@line[$n] eq 'TTC') {@line[$n] = 2;}
        if (@line[$n] eq 'TTA') {@line[$n] = 3;}
        if (@line[$n] eq 'TTG') {@line[$n] = 4;}
        if (@line[$n] eq 'CTT') {@line[$n] = 5;}
        if (@line[$n] eq 'CTC') {@line[$n] = 6;}
        if (@line[$n] eq 'CTA') {@line[$n] = 7;}
        if (@line[$n] eq 'CTG') {@line[$n] = 8;}
        if (@line[$n] eq 'ATT') {@line[$n] = 9;}
        if (@line[$n] eq 'ATC') {@line[$n] = 10;}
        if (@line[$n] eq 'ATA') {@line[$n] = 11;}
        if (@line[$n] eq 'ATG') {@line[$n] = 12;}
        if (@line[$n] eq 'GTT') {@line[$n] = 13;}
        if (@line[$n] eq 'GTC') {@line[$n] = 14;}
        if (@line[$n] eq 'GTA') {@line[$n] = 15;}
        if (@line[$n] eq 'GTG') {@line[$n] = 16;}
        if (@line[$n] eq 'TCT') {@line[$n] = 17;}
        if (@line[$n] eq 'TCC') {@line[$n] = 18;}
        if (@line[$n] eq 'TCA') {@line[$n] = 19;}
        if (@line[$n] eq 'TCG') {@line[$n] = 20;}
        if (@line[$n] eq 'CCT') {@line[$n] = 21;}
        if (@line[$n] eq 'CCC') {@line[$n] = 22;}
        if (@line[$n] eq 'CCA') {@line[$n] = 23;}
        if (@line[$n] eq 'CCG') {@line[$n] = 24;}
        if (@line[$n] eq 'ACT') {@line[$n] = 25;}
        if (@line[$n] eq 'ACC') {@line[$n] = 26;}
        if (@line[$n] eq 'ACA') {@line[$n] = 27;}
        if (@line[$n] eq 'ACG') {@line[$n] = 28;}
        if (@line[$n] eq 'GCT') {@line[$n] = 29;}
        if (@line[$n] eq 'GCC') {@line[$n] = 30;}
        if (@line[$n] eq 'GCA') {@line[$n] = 31;}
        if (@line[$n] eq 'GCG') {@line[$n] = 32;}
        if (@line[$n] eq 'TAT') {@line[$n] = 33;}
        if (@line[$n] eq 'TAC') {@line[$n] = 34;}
        if (@line[$n] eq 'TAA') {@line[$n] = 35;}
        if (@line[$n] eq 'TAG') {@line[$n] = 36;}
        if (@line[$n] eq 'CAT') {@line[$n] = 37;}
        if (@line[$n] eq 'CAC') {@line[$n] = 38;}
        if (@line[$n] eq 'CAA') {@line[$n] = 39;}
        if (@line[$n] eq 'CAG') {@line[$n] = 40;}
        if (@line[$n] eq 'AAT') {@line[$n] = 41;}
        if (@line[$n] eq 'AAC') {@line[$n] = 42;}
        if (@line[$n] eq 'AAA') {@line[$n] = 43;}
        if (@line[$n] eq 'AAG') {@line[$n] = 44;}
        if (@line[$n] eq 'GAT') {@line[$n] = 45;}
        if (@line[$n] eq 'GAC') {@line[$n] = 46;}
        if (@line[$n] eq 'GAA') {@line[$n] = 47;}
        if (@line[$n] eq 'GAG') {@line[$n] = 48;}
        if (@line[$n] eq 'TGT') {@line[$n] = 49;}
        if (@line[$n] eq 'TGC') {@line[$n] = 50;}
        if (@line[$n] eq 'TGA') {@line[$n] = 51;}
        if (@line[$n] eq 'TGG') {@line[$n] = 52;}
        if (@line[$n] eq 'CGT') {@line[$n] = 53;}
        if (@line[$n] eq 'CGC') {@line[$n] = 54;}
        if (@line[$n] eq 'CGA') {@line[$n] = 55;}
        if (@line[$n] eq 'CGG') {@line[$n] = 56;}
        if (@line[$n] eq 'AGT') {@line[$n] = 57;}
        if (@line[$n] eq 'AGC') {@line[$n] = 58;}
        if (@line[$n] eq 'AGA') {@line[$n] = 59;}
        if (@line[$n] eq 'AGG') {@line[$n] = 60;}
        if (@line[$n] eq 'GGT') {@line[$n] = 61;}
        if (@line[$n] eq 'GGC') {@line[$n] = 62;}
        if (@line[$n] eq 'GGA') {@line[$n] = 63;}
        if (@line[$n] eq 'GGG') {@line[$n] = 64;}
        
        chomp (@line[$n]);
        print (OUT1 "@line[$n]\t");
    }
    print (OUT1 "\n");
}
close (OUT1);
