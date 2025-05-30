open (IN1,"data.ffn");
@str1 = <IN1>;
$count=1;
$remove_polyc=0;
$remove_polys=0;
$total_polyc_m=0;
$total_polyc_s=0;
$codnum=0;

#number of species to which the two collapse sequence will be made
$collapse_group_num = $ARGV[0];

#for each species, get the order of the first and last sequence in the original alignment file
$x=0;
for ($n=0;$n<$collapse_group_num;$n++) {
    ${'st_group'.$n} = $ARGV[2*$n+1];
    ${'ed_group'.$n} = $ARGV[2*$n+2];
    $x = $x+2;
}

#number of species to which a single genome will be used for the analysis
$single_seq_num = $ARGV[$x+1];

#for each species, get the order of the single genome in the original alignment file
for ($n=0;$n<$single_seq_num;$n++) {
    ${'single_seq'.$n} = $ARGV[$x+2+$n];
}

#separate lines showing sequence name and nucleotide sequence
foreach $_ (@str1) {
    if ($_ =~/\>.*\n/) {
        $index[$count] = $_;
    }
    
    if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C||^N||^\./)) {
        $seq[$count] = $_;
        $count++;
    }
}

#convert nucleotide sequence string to array
for ($n=1;$n<=$count-1;$n++) {
    @{'SEQ'.$n} = split //, $seq[$n];
}

$x=0;

while ($x<$#SEQ1) {

    #for each of the species and each site of the sequences
    for ($n1=0;$n1<$collapse_group_num;$n1++) {

        $collapse_index_A = 2*$n1+1;
        $collapse_index_B = 2*$n1+2;
        
        #get the nucleotide state of the first sequence of the species and put it to $refpos1 and $refpos2
        $polyc=1;
        $refpos1 = @{'SEQ'.${'st_group'.$n1}}[$x];
        $refpos2 = @{'SEQ'.${'st_group'.$n1}}[$x];

        #check the second to the last sequence of the species
        $polys=0;
        for ($n2=${'st_group'.$n1}+1;$n2<=${'ed_group'.$n1};$n2++) {
            $pos = @{'SEQ'.$n2}[$x];
            
            #check the nucleotide state is different from both of $refpos1 and $refpos2, increment $polys value and put the nucleotide sate to $refpos2
            if (($refpos1 ne $pos) && ($refpos2 ne $pos)) {
                $polys++;
                $refpos2 = $pos;
            }
        }
        #site with >2 states should be filtered beforehand, so if there is, output the error message
        if ($polys >= 2) {
            printf "ERROR: more than two states at $x th polymorphic site of $n1 th species\n";
        }
        
        #$polys == 1 menas there are 2 states at the site
        if ($polys == 1) {
            
            $r = int (rand(10));

            #one of the collapse sequence has the first state and the other has the second state
            #which collapse sequence has which sate is decided randomly
            #'nSEQ' are the output sequences
            if ($r<=4) {
                @{'nSEQ'.$collapse_index_A}[$x] = $refpos1;
                @{'nSEQ'.$collapse_index_B}[$x] = $refpos2;
            }
            elsif ($r>4) {
                @{'nSEQ'.$collapse_index_A}[$x] = $refpos2;
                @{'nSEQ'.$collapse_index_B}[$x] = $refpos1;
            }
        }
    
        #$polys == 0 menas the site is monomorphic
        if ($polys == 0) {
            
            #assign the monomorphic state to both of the collapse sequences
            @{'nSEQ'.$collapse_index_A}[$x] = $refpos1;
            @{'nSEQ'.$collapse_index_B}[$x] = $refpos1;
        }
    }
    
    #for the species with a single genome, just copy the input sequence to the output 'nSEQ'
    $start = 2*$collapse_group_num+1;
    
    for ($n3=0;$n3<$single_seq_num;$n3++) {
        $ref= $start+$n3;
        @{'nSEQ'.$ref}[$x] = @{'SEQ'.${'single_seq'.$n3}}[$x];
    }

    
    $x++;
}

# output the 'nSEQ' to the file
open (OUT1, ">>ndata.ffn");

$start = 2*$collapse_group_num+1;
for ($n3=0;$n3<$single_seq_num;$n3++) {
    
    print (OUT1 "\>single\_seq\_$n3\n");
    $ref = $start+$n3;
    $x=0;
    while ($x<=$#nSEQ1) {
        if (@{'nSEQ'.$ref} !~ /\s/) {
            print (OUT1 "@{'nSEQ'.$ref}[$x]");
        }
        $x++;
    }
    print (OUT1 "\n");
}


for ($n3=0;$n3<$collapse_group_num;$n3++) {
    
    print (OUT1 "\>collapse\_seq\_$n3\_A\n");
    $refA = 2*$n3+1;
    $x=0;
    while ($x<=$#nSEQ1) {
        if (@{'nSEQ'.$refA} !~ /\s/) {
            print (OUT1 "@{'nSEQ'.$refA}[$x]");
        }
        $x++;
    }
    print (OUT1 "\n");
    print (OUT1 "\>collapse\_seq\_$n3\_B\n");
    $refB = 2*$n3+2;
    $x=0;
    while ($x<=$#nSEQ1) {
        if (@{'nSEQ'.$refB} !~ /\s/) {
            print (OUT1 "@{'nSEQ'.$refB}[$x]");
        }
        $x++;
    }
    print (OUT1 "\n");
}

close (OUT1);


