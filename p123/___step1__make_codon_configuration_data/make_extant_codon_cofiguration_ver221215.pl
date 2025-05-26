open (IN1,"data/seq_dat_pos1.txt");
@str1 = <IN1>;
open (IN2,"data/seq_dat_pos2.txt");
@str2 = <IN2>;
open (IN3,"data/seq_dat_pos3.txt");
@str3 = <IN3>;

$codnum=0;

$count_pos1=0;
foreach $_ (@str1) {
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C||^N||^\./)) {
        $seq_pos1[$count_pos1] = $_;
        $count_pos1++;
    }
}
$count_pos2=0;
foreach $_ (@str2) {
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C||^N||^\./)) {
        $seq_pos2[$count_pos2] = $_;
        $count_pos2++;
    }
}
$count_pos3=0;
foreach $_ (@str3) {
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C||^N||^\./)) {
        $seq_pos3[$count_pos3] = $_;
        $count_pos3++;
    }
}

if (($count_pos1 != $count_pos2)||($count_pos2 != $count_pos3)||($count_pos3 != $count_pos1)) {
    printf "ERROR: seq length of three positions are not equal\n";
}

for ($n=0;$n<=$count_pos1-1;$n++) {
    @{'SEQ_pos1'.$n} = split //, $seq_pos1[$n];
    @{'SEQ_pos2'.$n} = split //, $seq_pos2[$n];
    @{'SEQ_pos3'.$n} = split //, $seq_pos3[$n];
}

$x=0;

while ($x<$#SEQ_pos10) {
    $check=1;
    
    $refpos1_mel = @{'SEQ_pos1'.0}[$x];
    $refpos2_mel = @{'SEQ_pos2'.0}[$x];
    $refpos3_mel = @{'SEQ_pos3'.0}[$x];
    $refcodon1_mel = "$refpos1_mel$refpos2_mel$refpos3_mel";
    $refcodon2_mel = "XXX";
    $freq_ref1_mel = 1;
    $freq_ref2_mel = 0;
    
    for ($n=1;$n<=13;$n++) {
        $pos1_mel = @{'SEQ_pos1'.$n}[$x];
        $pos2_mel = @{'SEQ_pos2'.$n}[$x];
        $pos3_mel = @{'SEQ_pos3'.$n}[$x];
        $codon_mel = "$pos1_mel$pos2_mel$pos3_mel";
        
        if ($codon_mel eq $refcodon1_mel) {
            $freq_ref1_mel++;
        }
        if ($codon_mel ne $refcodon1_mel) {
            if ($refcodon2_mel eq "XXX") {
                $refcodon2_mel = $codon_mel;
                $freq_ref2_mel++;
            }
            elsif ($refcodon2_mel ne "XXX") {
                if ($refcodon2_mel eq $codon_mel) {
                    $freq_ref2_mel++;
                }
                elsif ($refcodon2_mel ne $codon_mel) {
                    printf "ERROR: There are more than two codons in Dmel starting from $x th site\n";
                }
            }
        }
    }
    $refpos1_sim = @{'SEQ_pos1'.14}[$x];
    $refpos2_sim = @{'SEQ_pos2'.14}[$x];
    $refpos3_sim = @{'SEQ_pos3'.14}[$x];
    $refcodon1_sim = "$refpos1_sim$refpos2_sim$refpos3_sim";
    $refcodon2_sim = "XXX";
    $freq_ref1_sim = 1;
    $freq_ref2_sim = 0;
    
    for ($n=15;$n<=34;$n++) {
        $pos1_sim = @{'SEQ_pos1'.$n}[$x];
        $pos2_sim = @{'SEQ_pos2'.$n}[$x];
        $pos3_sim = @{'SEQ_pos3'.$n}[$x];
        $codon_sim = "$pos1_sim$pos2_sim$pos3_sim";
        
        if ($codon_sim eq $refcodon1_sim) {
            $freq_ref1_sim++;
        }
        if ($codon_sim ne $refcodon1_sim) {
            if ($refcodon2_sim eq "XXX") {
                $refcodon2_sim = $codon_sim;
                $freq_ref2_sim++;
            }
            elsif ($refcodon2_sim ne "XXX") {
                if ($refcodon2_sim eq $codon_sim) {
                    $freq_ref2_sim++;
                }
                elsif ($refcodon2_sim ne $codon_sim) {
                    printf "ERROR: There are more than two codons in Dsim starting from $x th site\n";
                }
            }
        }
    }
    
    $refpos1_yak = @{'SEQ_pos1'.35}[$x];
    $refpos2_yak = @{'SEQ_pos2'.35}[$x];
    $refpos3_yak = @{'SEQ_pos3'.35}[$x];
    $refcodon1_yak = "$refpos1_yak$refpos2_yak$refpos3_yak";
    $refcodon2_yak = "XXX";
    $freq_ref1_yak = 1;
    $freq_ref2_yak = 0;
    
    for ($n=36;$n<=35;$n++) {
        $pos1_yak = @{'SEQ_pos1'.$n}[$x];
        $pos2_yak = @{'SEQ_pos2'.$n}[$x];
        $pos3_yak = @{'SEQ_pos3'.$n}[$x];
        $codon_yak = "$pos1_yak$pos2_yak$pos3_yak";
        
        if ($codon_yak eq $refcodon1_yak) {
            $freq_ref1_yak++;
        }
        if ($codon_yak ne $refcodon1_yak) {
            if ($refcodon2_yak eq "XXX") {
                $refcodon2_yak = $codon_yak;
                $freq_ref2_yak++;
            }
            if ($refcodon2_yak ne "XXX") {
                if ($refcodon2_yak eq $codon_yak) {
                    $freq_ref2_yak++;
                }
                if ($refcodon2_yak ne $codon_yak) {
                    printf "ERROR: There are more than two codons in Dyak starting from $x th site\n";
                }
            }
        }
    }
    
    $refpos1_ere = @{'SEQ_pos1'.36}[$x];
    $refpos2_ere = @{'SEQ_pos2'.36}[$x];
    $refpos3_ere = @{'SEQ_pos3'.36}[$x];
    $refcodon1_ere = "$refpos1_ere$refpos2_ere$refpos3_ere";
    $refcodon2_ere = "XXX";
    $freq_ref1_ere = 1;
    $freq_ref2_ere = 0;
    
    for ($n=37;$n<=36;$n++) {
        $pos1_ere = @{'SEQ_pos1'.$n}[$x];
        $pos2_ere = @{'SEQ_pos2'.$n}[$x];
        $pos3_ere = @{'SEQ_pos3'.$n}[$x];
        $codon_ere = "$pos1_ere$pos2_ere$pos3_ere";
        
        if ($codon_ere eq $refcodon1_ere) {
            $freq_ref1_ere++;
        }
        if ($codon_ere ne $refcodon1_ere) {
            if ($refcodon2_ere eq "XXX") {
                $refcodon2_ere = $codon_ere;
                $freq_ref2_ere++;
            }
            if ($refcodon2_ere ne "XXX") {
                if ($refcodon2_ere eq $codon_ere) {
                    $freq_ref2_ere++;
                }
                if ($refcodon2_ere ne $codon_ere) {
                    printf "ERROR: There are more than two codons in Dere starting from $x th site\n";
                }
            }
        }
    }
    
    
    open (OUT1, ">>codon_configuration_at_extant_nodes.txt");
    
    print (OUT1 "$x\t");
    
    print (OUT1 "$refcodon1_mel\t$freq_ref1_mel\t");
    if ($refcodon2_mel ne "XXX") {
        print (OUT1 "$refcodon2_mel\t$freq_ref2_mel\t");
    }
    if ($refcodon2_mel eq "XXX") {
        print (OUT1 "none\t0\t");
    }
    print (OUT1 "$refcodon1_sim\t$freq_ref1_sim\t");
    if ($refcodon2_sim ne "XXX") {
        print (OUT1 "$refcodon2_sim\t$freq_ref2_sim\t");
    }
    if ($refcodon2_sim eq "XXX") {
        print (OUT1 "none\t0\t");
    }
    print (OUT1 "$refcodon1_yak\t$freq_ref1_yak\t");
    if ($refcodon2_yak ne "XXX") {
        print (OUT1 "$refcodon2_yak\t$freq_ref2_yak\t");
    }
    if ($refcodon2_yak eq "XXX") {
        print (OUT1 "none\t0\t");
    }
    print (OUT1 "$refcodon1_ere\t$freq_ref1_ere\t");
    if ($refcodon2_ere ne "XXX") {
        print (OUT1 "$refcodon2_ere\t$freq_ref2_ere\t");
    }
    if ($refcodon2_ere eq "XXX") {
        print (OUT1 "none\t0\t");
    }
    print (OUT1 "\n");
    
    $x++;
}
 
