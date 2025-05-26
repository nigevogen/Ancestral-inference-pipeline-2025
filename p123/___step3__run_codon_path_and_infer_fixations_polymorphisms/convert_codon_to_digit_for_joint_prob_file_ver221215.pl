open (IN1,"joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons.txt");
@str1 = <IN1>;

open (OUT1, ">>joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons_codon_digit.txt");

foreach $_ (@str1) {
    @line = split /\t/, $_;
    $codon1 = @line[1];
    $codon2 = @line[2];
    $codon3 = @line[3];
    $codon4 = @line[4];
    
    for ($n=1;$n<=4;$n++) {
        if (${'codon'.$n} eq 'TTT') {${'codon'.$n} = 1;}
        if (${'codon'.$n} eq 'TTC') {${'codon'.$n} = 2;}
        if (${'codon'.$n} eq 'TTA') {${'codon'.$n} = 3;}
        if (${'codon'.$n} eq 'TTG') {${'codon'.$n} = 4;}
        if (${'codon'.$n} eq 'CTT') {${'codon'.$n} = 5;}
        if (${'codon'.$n} eq 'CTC') {${'codon'.$n} = 6;}
        if (${'codon'.$n} eq 'CTA') {${'codon'.$n} = 7;}
        if (${'codon'.$n} eq 'CTG') {${'codon'.$n} = 8;}
        if (${'codon'.$n} eq 'ATT') {${'codon'.$n} = 9;}
        if (${'codon'.$n} eq 'ATC') {${'codon'.$n} = 10;}
        if (${'codon'.$n} eq 'ATA') {${'codon'.$n} = 11;}
        if (${'codon'.$n} eq 'ATG') {${'codon'.$n} = 12;}
        if (${'codon'.$n} eq 'GTT') {${'codon'.$n} = 13;}
        if (${'codon'.$n} eq 'GTC') {${'codon'.$n} = 14;}
        if (${'codon'.$n} eq 'GTA') {${'codon'.$n} = 15;}
        if (${'codon'.$n} eq 'GTG') {${'codon'.$n} = 16;}
        if (${'codon'.$n} eq 'TCT') {${'codon'.$n} = 17;}
        if (${'codon'.$n} eq 'TCC') {${'codon'.$n} = 18;}
        if (${'codon'.$n} eq 'TCA') {${'codon'.$n} = 19;}
        if (${'codon'.$n} eq 'TCG') {${'codon'.$n} = 20;}
        if (${'codon'.$n} eq 'CCT') {${'codon'.$n} = 21;}
        if (${'codon'.$n} eq 'CCC') {${'codon'.$n} = 22;}
        if (${'codon'.$n} eq 'CCA') {${'codon'.$n} = 23;}
        if (${'codon'.$n} eq 'CCG') {${'codon'.$n} = 24;}
        if (${'codon'.$n} eq 'ACT') {${'codon'.$n} = 25;}
        if (${'codon'.$n} eq 'ACC') {${'codon'.$n} = 26;}
        if (${'codon'.$n} eq 'ACA') {${'codon'.$n} = 27;}
        if (${'codon'.$n} eq 'ACG') {${'codon'.$n} = 28;}
        if (${'codon'.$n} eq 'GCT') {${'codon'.$n} = 29;}
        if (${'codon'.$n} eq 'GCC') {${'codon'.$n} = 30;}
        if (${'codon'.$n} eq 'GCA') {${'codon'.$n} = 31;}
        if (${'codon'.$n} eq 'GCG') {${'codon'.$n} = 32;}
        if (${'codon'.$n} eq 'TAT') {${'codon'.$n} = 33;}
        if (${'codon'.$n} eq 'TAC') {${'codon'.$n} = 34;}
        if (${'codon'.$n} eq 'TAA') {${'codon'.$n} = 35;}
        if (${'codon'.$n} eq 'TAG') {${'codon'.$n} = 36;}
        if (${'codon'.$n} eq 'CAT') {${'codon'.$n} = 37;}
        if (${'codon'.$n} eq 'CAC') {${'codon'.$n} = 38;}
        if (${'codon'.$n} eq 'CAA') {${'codon'.$n} = 39;}
        if (${'codon'.$n} eq 'CAG') {${'codon'.$n} = 40;}
        if (${'codon'.$n} eq 'AAT') {${'codon'.$n} = 41;}
        if (${'codon'.$n} eq 'AAC') {${'codon'.$n} = 42;}
        if (${'codon'.$n} eq 'AAA') {${'codon'.$n} = 43;}
        if (${'codon'.$n} eq 'AAG') {${'codon'.$n} = 44;}
        if (${'codon'.$n} eq 'GAT') {${'codon'.$n} = 45;}
        if (${'codon'.$n} eq 'GAC') {${'codon'.$n} = 46;}
        if (${'codon'.$n} eq 'GAA') {${'codon'.$n} = 47;}
        if (${'codon'.$n} eq 'GAG') {${'codon'.$n} = 48;}
        if (${'codon'.$n} eq 'TGT') {${'codon'.$n} = 49;}
        if (${'codon'.$n} eq 'TGC') {${'codon'.$n} = 50;}
        if (${'codon'.$n} eq 'TGA') {${'codon'.$n} = 51;}
        if (${'codon'.$n} eq 'TGG') {${'codon'.$n} = 52;}
        if (${'codon'.$n} eq 'CGT') {${'codon'.$n} = 53;}
        if (${'codon'.$n} eq 'CGC') {${'codon'.$n} = 54;}
        if (${'codon'.$n} eq 'CGA') {${'codon'.$n} = 55;}
        if (${'codon'.$n} eq 'CGG') {${'codon'.$n} = 56;}
        if (${'codon'.$n} eq 'AGT') {${'codon'.$n} = 57;}
        if (${'codon'.$n} eq 'AGC') {${'codon'.$n} = 58;}
        if (${'codon'.$n} eq 'AGA') {${'codon'.$n} = 59;}
        if (${'codon'.$n} eq 'AGG') {${'codon'.$n} = 60;}
        if (${'codon'.$n} eq 'GGT') {${'codon'.$n} = 61;}
        if (${'codon'.$n} eq 'GGC') {${'codon'.$n} = 62;}
        if (${'codon'.$n} eq 'GGA') {${'codon'.$n} = 63;}
        if (${'codon'.$n} eq 'GGG') {${'codon'.$n} = 64;}
    }
    

    print (OUT1 "@line[0]\t$codon1\t$codon2\t$codon3\t$codon4\t@line[5]");

}
close (OUT1);
