#The purpose of this code is to weight the estimated ancestral probabilities using the estimated SFS of the previous round of iterative BTW

open (IN1,"collapse_seq.txt"); #collapse_alignment which was the input of BASEML
@str1 = <IN1>;
open (IN2,"fltrst_re"); #output of BASEML which shows the estimated ancestral probabilities and tree information
@str2 = <IN2>;
open (IN3,"original_aln.txt"); #original_alignment which was the input to make collapse sequences
@str3 = <IN3>;
open (IN5,"target_node"); #file showing the position of the target node (ancestor of the two collapse sequences of the focused species (Xth species in the shell script)
@str5 = <IN5>;

open (IN4,"estimated_frequency_spectrum.txt"); #estimated SFS of the previous round
@str4 = <IN4>;

# read input parameters

$total_seq_num = $ARGV[0]; #number of the sequences in the collapse_alignment
$num_internal_nodes = $ARGV[1]; #number of internal nodes in the input tree of BASEML
$sample_start = $ARGV[3]; #order of the first sequence of the focused species in the original_alignemnt
$sample_end = $ARGV[4]; #order of the last sequence of the focused species in the original_alignemnt
$collapse_start = $ARGV[17]; #order of the first collapse sequence of the focused species in the collapse_alignemnt
$collapse_end = $ARGV[18]; #order of the second collapse sequence of the focused species in the collapse_alignemnt
$target_node = @str5[0];
chomp ($target_node);

# calculate the number of the sequences of the focused species in original_alignment
$sample_num = $sample_end - $sample_start + 1;

$mutation_category_num = $ARGV[2]; #number of mutation category having different expected SFS (SFS category)

#read SFS of each of the 12 mutation category
@TC = split /\t/, @str4[0];
@TA = split /\t/, @str4[1];
@TG = split /\t/, @str4[2];
@CT = split /\t/, @str4[3];
@CA = split /\t/, @str4[4];
@CG = split /\t/, @str4[5];
@AT = split /\t/, @str4[6];
@AC = split /\t/, @str4[7];
@AG = split /\t/, @str4[8];
@GT = split /\t/, @str4[9];
@GC = split /\t/, @str4[10];
@GA = split /\t/, @str4[11];

$num_TC=0;
$num_TA=0;
$num_TG=0;
$num_CT=0;
$num_CA=0;
$num_CG=0;
$num_AT=0;
$num_AC=0;
$num_AG=0;
$num_GT=0;
$num_GC=0;
$num_GA=0;

# sum the numbers of polymorphic mutations in one mutation category
for ($n1=1;$n1<=$mutation_category_num;$n1++) {
    for ($n2=1;$n2<=$sample_num;$n2++) {
        @{'num_category'.$n1}[$n2] = 0;
    }
    
    # check the SFS category of each of the 12 mutation categories from the input parameters
    $cate_check_TC = $ARGV[5];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[5]}[$n2] = @{'num_category'.$ARGV[5]}[$n2]+@TC[$n2];
    }
    $cate_check_TA = $ARGV[6];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[6]}[$n2] = @{'num_category'.$ARGV[6]}[$n2]+@TA[$n2];
    }
    $cate_check_TG = $ARGV[7];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[7]}[$n2] = @{'num_category'.$ARGV[7]}[$n2]+@TG[$n2];
    }
    $cate_check_CT = $ARGV[8];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[8]}[$n2] = @{'num_category'.$ARGV[8]}[$n2]+@CT[$n2];
    }
    $cate_check_CA = $ARGV[9];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[9]}[$n2] = @{'num_category'.$ARGV[9]}[$n2]+@CA[$n2];
    }
    $cate_check_CG = $ARGV[10];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[10]}[$n2] = @{'num_category'.$ARGV[10]}[$n2]+@CG[$n2];
    }
    $cate_check_AT = $ARGV[11];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[11]}[$n2] = @{'num_category'.$ARGV[11]}[$n2]+@AT[$n2];
    }
    $cate_check_AC = $ARGV[12];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[12]}[$n2] = @{'num_category'.$ARGV[12]}[$n2]+@AC[$n2];
    }
    $cate_check_AG = $ARGV[13];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[13]}[$n2] = @{'num_category'.$ARGV[13]}[$n2]+@AG[$n2];
    }
    $cate_check_GT = $ARGV[14];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[14]}[$n2] = @{'num_category'.$ARGV[14]}[$n2]+@GT[$n2];
    }
    $cate_check_GC = $ARGV[15];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[15]}[$n2] = @{'num_category'.$ARGV[15]}[$n2]+@GC[$n2];
    }
    $cate_check_GA = $ARGV[16];
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        @{'num_category'.$ARGV[16]}[$n2] = @{'num_category'.$ARGV[16]}[$n2]+@GA[$n2];
    }
}

#output the SFS category
printf "TC is in category $cate_check_TC\n";
printf "TA is in category $cate_check_TA\n";
printf "TG is in category $cate_check_TG\n";
printf "CT is in category $cate_check_CT\n";
printf "CA is in category $cate_check_CA\n";
printf "CG is in category $cate_check_CG\n";
printf "AT is in category $cate_check_AT\n";
printf "AC is in category $cate_check_AC\n";
printf "AG is in category $cate_check_AG\n";
printf "GT is in category $cate_check_GT\n";
printf "GC is in category $cate_check_GC\n";
printf "GA is in category $cate_check_GA\n";

# calculate the SFSest of each mutation category
# Be careful for the order of each frequency class in SFS below, which is (n, n-1, n-2, n-3, ..., 1, 0)
for ($n1=1;$n1<=$mutation_category_num;$n1++) {
    $total = 0;
    for ($n2=1;$n2<=$sample_num-1;$n2++) {
        $total = $total + @{'num_category'.$n1}[$n2];
    }
    if ($total > 0) {
        $m=1;
        for ($n2=$sample_num-1;$n2>=1;$n2--) {
            @{'freq_category'.$n1}[$m] = @{'num_category'.$n1}[$n2]/$total;
            $m++;
        }
    }
    if ($total == 0) {
        for ($n2=1;$n2<=$sample_num-1;$n2++) {
            @{'freq_category'.$n1}[$n2] = 0;
        }
    }
    # put 1 and 0 for frequency 1 and 0
    @{'freq_category'.$n1}[0] = 0;
    @{'freq_category'.$n1}[$sample_num] = 1;
}

#output the estimated SFS of each SFS category
for ($n1=1;$n1<=$mutation_category_num;$n1++) {
    printf "@{'freq_category'.$n1}\n";
}

# make array to show each site of each sequences in the collapse_alignment
# first put a sequence into string and then convert to array using split function
# we can access each site of each sequence as @{'SEQ'.$m}[$n] (nth site of mth sequence)
for ($n1=1;$n1<=$total_seq_num;$n1++) {
    $n2 = $n1-1;
    ${'SEQ'.$n1} = @str1[2*$n2+1];
    @{'SEQ'.$n1} = split //, ${'SEQ'.$n1};
}

# make array to show each site of each sequences of the focused species in the original_alignment
$n3 = $sample_start-1;
for ($n1=1;$n1<=$sample_num;$n1++) {
    ${'mSEQ'.$n1} = @str3[2*$n3+1];
    @{'mSEQ'.$n1} = split //, ${'mSEQ'.$n1};
    $n3++;
}

# output file
open (OUT1, ">>anc_site_probs.txt");
open (OUT2, ">>anc_site_probs_all_node.txt");

$n1=0;

# for each site (from 1st to the end)
for ($n=0;$n<=$#SEQ1-1;$n++) {
    $site_pattern = "";
    # generate the nucleotide states configuration of the sequences in the collapse_alignment
    for ($x=1;$x<=$total_seq_num;$x++) {
        $site_pattern = $site_pattern.@{'SEQ'.$x}[$n];
    }
    
    # check the polymorphism state of the focused species at the current site
    # there should be 2 polymoprhic states (no more than 2), # of the 1st state is $nref1 and # of the 2nd state is $nref2
    $nref1=0;
    $nref2=0;
    
    # the 1st state is that in the 1st sequence of mSEQ = @{'mSEQ'.1}
    $site = $n;
    $ref1 = @mSEQ1[$site];
    # first we assign 'nan' for the 2nd state
    $ref2 = 'nan';
    # the 1st sequence of mSEQ has the 1st state, so # of the 1st state starts from 1
    $nref1 = 1;
    #check the 2nd to the last sequences and update the count of the 1st and 2nd states
    for ($n2=2;$n2<=$sample_num;$n2++) {
        
        if (@{'mSEQ'.$n2}[$site] eq $ref1) {
            $nref1++;
        }
        if (@{'mSEQ'.$n2}[$site] ne $ref1) {
            $ref2 = @{'mSEQ'.$n2}[$site];
            $nref2++;
        }
        #output error if there are more than two states
        if ((@{'mSEQ'.$n2}[$site] ne $ref1) && (@{'mSEQ'.$n2}[$site] ne $ref2)) {
            printf "ERROR: there are more than two states at $n th site: $ref1 $ref2 @{'mSEQ'.$n2}[$site]\n";
        }
    }
    $freq_ref1 = $nref1;
    $freq_ref2 = $nref2;
    
    $n1++;
    
    #read fltrst_re and find the line showing the ancestral states configuration of the current states configuration of the current site
    #for the nth site of the alignment, nth line of the fltrst shows this info
    $line = @str2[$site];
    if ($line =~ /^(\d+)\s+1\s+$site_pattern:\s(\w+)/) {
        
        print (OUT1 "$site\t");
        #output the info of $site, frequencies of 1st and 2nd polymorphic states (freq of the state that the 1st collapse seq has comes first)
        if (@{'SEQ'.$collapse_start}[$site] eq $ref1) {
            print (OUT2 "$site\t$freq_ref1\t$freq_ref2\t$site_pattern\:\t\t");
        }
        elsif (@{'SEQ'.$collapse_start}[$site] eq $ref2) {
            print (OUT2 "$site\t$freq_ref2\t$freq_ref1\t$site_pattern\:\t\t");
        }
        
        $num = $1;
        # find the inferred state at the ancestral node of the two collapse sequences of the focused species
        $preref = $2;
        @preref = split //, $preref;
        $ref_m = @preref[$target_node-1];
        
        # check the frequency of the two polymorphic states (frequency of the state consistent with the ancestral state in the first inference is stored as $freq)
        if ($ref_m eq $ref1) {
            $freq = $freq_ref1;
        }
        if ($ref_m eq $ref2) {
            $freq = $freq_ref2;
        }
        
        $total = 0;
        # weight the probability using the expected SFS (1ST STEP: calculate the total of the weighted probabilities)
        # for each of the inferred ancestral configuration...
        while ($line =~ /(\w+)\s+\((\d\.\d+)\)/g) {
            $preref = $1;
            
            # the last component matching to the above pattern is the total probability, so skip this matching
            if ($preref ne "total") {
                @preref = split //, $preref;
                $ref_m = @preref[$target_node-1];
                
                # if the ancestral state is consistent with the state $ref1, the expected proportion of the derived mutation in SFSne is @SFSest[$freq_ref1]. Again please be careful for the order of the SFS defined at the beginning
                if ($ref_m eq $ref1) {
                    
                    # weight for each ancestral and derived nucleotide pair
                    if ($ref_m eq "A") {
                        # $ref2 eq "nan" means monomorphic, so 1 is multiplied here
                        if ($ref2 eq "nan") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_AT}[$freq_ref1]);
                        }
                        if ($ref2 eq "T") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_AT}[$freq_ref1]);
                        }
                        if ($ref2 eq "G") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_AG}[$freq_ref1]);
                        }
                        if ($ref2 eq "C") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_AC}[$freq_ref1]);
                        }
                    }
                    if ($ref_m eq "T") {
                        if ($ref2 eq "A") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_TA}[$freq_ref1]);
                        }
                        if ($ref2 eq "nan") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_TA}[$freq_ref1]);
                        }
                        if ($ref2 eq "G") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_TG}[$freq_ref1]);
                        }
                        if ($ref2 eq "C") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_TC}[$freq_ref1]);
                        }
                    }
                    if ($ref_m eq "G") {
                        if ($ref2 eq "A") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_GA}[$freq_ref1]);
                        }
                        if ($ref2 eq "T") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_GT}[$freq_ref1]);
                        }
                        if ($ref2 eq "nan") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_GC}[$freq_ref1]);
                        }
                        if ($ref2 eq "C") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_GC}[$freq_ref1]);
                        }
                    }
                    if ($ref_m eq "C") {
                        if ($ref2 eq "A") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_CA}[$freq_ref1]);
                        }
                        if ($ref2 eq "T") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_CT}[$freq_ref1]);
                        }
                        if ($ref2 eq "G") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_CG}[$freq_ref1]);
                        }
                        if ($ref2 eq "nan") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_CG}[$freq_ref1]);
                        }
                    }
                }
                
                # if the ancestral state is consistent with the state $ref2, the expected proportion of the derived mutation in SFSne is @SFSest[$freq_ref2].
                if ($ref_m eq $ref2) {
                    if ($ref_m eq "A") {
                        if ($ref1 eq "A") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_AT}[$freq_ref2]);
                        }
                        if ($ref1 eq "T") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_AT}[$freq_ref2]);
                        }
                        if ($ref1 eq "G") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_AG}[$freq_ref2]);
                        }
                        if ($ref1 eq "C") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_AC}[$freq_ref2]);
                        }
                    }
                    if ($ref_m eq "T") {
                        if ($ref1 eq "A") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_TA}[$freq_ref2]);
                        }
                        if ($ref1 eq "T") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_TA}[$freq_ref2]);
                        }
                        if ($ref1 eq "G") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_TG}[$freq_ref2]);
                        }
                        if ($ref1 eq "C") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_TC}[$freq_ref2]);
                        }
                    }
                    if ($ref_m eq "G") {
                        if ($ref1 eq "A") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_GA}[$freq_ref2]);
                        }
                        if ($ref1 eq "T") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_GT}[$freq_ref2]);
                        }
                        if ($ref1 eq "G") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_GC}[$freq_ref2]);
                        }
                        if ($ref1 eq "C") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_GC}[$freq_ref2]);
                        }
                    }
                    if ($ref_m eq "C") {
                        if ($ref1 eq "A") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_CA}[$freq_ref2]);
                        }
                        if ($ref1 eq "T") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_CT}[$freq_ref2]);
                        }
                        if ($ref1 eq "G") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_CG}[$freq_ref2]);
                        }
                        if ($ref1 eq "C") {
                            $total = $total + ($2*@{'freq_category'.$cate_check_CG}[$freq_ref2]);
                        }
                    }
                }
                
                # if the ancestral state does not exist in the original_alignment
                if (($ref_m ne $ref1) && ($ref_m ne $ref2)) {
                    
                    # put probability 0 if the ancestral state is different from both of the polymorphic states
                    if (($ref1 ne "nan") && ($ref2 ne "nan")) {
                        $total = $total + 0.00;
                    }
                    # put unweighted probability if the site is monomorphic
                    if (($ref1 eq "nan") || ($ref2 eq "nan")) {
                        $total = $total + 0.00;
                    }
                }
            }
        }
        
        
        # weight the probability using the expected SFS (2ND: calculate the weighted probabilities for each ancestral configuration)
        # basically the repetetion of the above process
        if ($total > 0) {
            while ($line =~ /(\w+)\s+\((\d\.\d+)\)/g) {
                $preref = $1;
                @preref = split //, $preref;
                $ref_m = @preref[$target_node-1];
                
                for ($n7=0;$n7<=$num_internal_nodes-1;$n7++) {
                    ${'site'.$n7} = @preref[$n7];
                }
                
                # weight for each ancestral and derived nucleotide pair. This time, the weighted probability is devided by $total calculated above to make the total probability to be 1.
                if ($ref_m eq $ref1) {
                    if ($ref_m eq "A") {
                        if ($ref2 eq "nan") {
                            $X = ($2*@{'freq_category'.$cate_check_AT}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "T") {
                            $X = ($2*@{'freq_category'.$cate_check_AT}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "G") {
                            $X = ($2*@{'freq_category'.$cate_check_AG}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "C") {
                            $X = ($2*@{'freq_category'.$cate_check_AC}[$freq_ref1])/$total;
                        }
                    }
                    if ($ref_m eq "T") {
                        if ($ref2 eq "A") {
                            $X = ($2*@{'freq_category'.$cate_check_TA}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "nan") {
                            $X = ($2*@{'freq_category'.$cate_check_TA}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "G") {
                            $X = ($2*@{'freq_category'.$cate_check_TG}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "C") {
                            $X = ($2*@{'freq_category'.$cate_check_TC}[$freq_ref1])/$total;
                        }
                    }
                    if ($ref_m eq "G") {
                        if ($ref2 eq "A") {
                            $X = ($2*@{'freq_category'.$cate_check_GA}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "T") {
                            $X = ($2*@{'freq_category'.$cate_check_GT}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "nan") {
                            $X = ($2*@{'freq_category'.$cate_check_GC}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "C") {
                            $X = ($2*@{'freq_category'.$cate_check_GC}[$freq_ref1])/$total;
                        }
                    }
                    if ($ref_m eq "C") {
                        if ($ref2 eq "A") {
                            $X = ($2*@{'freq_category'.$cate_check_CA}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "T") {
                            $X = ($2*@{'freq_category'.$cate_check_CT}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "G") {
                            $X = ($2*@{'freq_category'.$cate_check_CG}[$freq_ref1])/$total;
                        }
                        if ($ref2 eq "nan") {
                            $X = ($2*@{'freq_category'.$cate_check_CG}[$freq_ref1])/$total;
                        }
                    }
                }
                
                if ($ref_m eq $ref2) {
                    if ($ref_m eq "A") {
                        if ($ref1 eq "A") {
                            $X = ($2*@{'freq_category'.$cate_check_AT}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "T") {
                            $X = ($2*@{'freq_category'.$cate_check_AT}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "G") {
                            $X = ($2*@{'freq_category'.$cate_check_AG}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "C") {
                            $X = ($2*@{'freq_category'.$cate_check_AC}[$freq_ref2])/$total;
                        }
                    }
                    if ($ref_m eq "T") {
                        if ($ref1 eq "A") {
                            $X = ($2*@{'freq_category'.$cate_check_TA}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "T") {
                            $X = ($2*@{'freq_category'.$cate_check_TA}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "G") {
                            $X = ($2*@{'freq_category'.$cate_check_TG}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "C") {
                            $X = ($2*@{'freq_category'.$cate_check_TC}[$freq_ref2])/$total;
                        }
                    }
                    if ($ref_m eq "G") {
                        if ($ref1 eq "A") {
                            $X = ($2*@{'freq_category'.$cate_check_GA}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "T") {
                            $X = ($2*@{'freq_category'.$cate_check_GT}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "G") {
                            $X = ($2*@{'freq_category'.$cate_check_GC}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "C") {
                            $X = ($2*@{'freq_category'.$cate_check_GC}[$freq_ref2])/$total;
                        }
                    }
                    if ($ref_m eq "C") {
                        if ($ref1 eq "A") {
                            $X = ($2*@{'freq_category'.$cate_check_CA}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "T") {
                            $X = ($2*@{'freq_category'.$cate_check_CT}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "G") {
                            $X = ($2*@{'freq_category'.$cate_check_CG}[$freq_ref2])/$total;
                        }
                        if ($ref1 eq "C") {
                            $X = ($2*@{'freq_category'.$cate_check_CG}[$freq_ref2])/$total;
                        }
                    }
                }
                
                if (($ref_m ne $ref1) && ($ref_m ne $ref2)) {
                    
                    if (($ref1 ne "nan") && ($ref2 ne "nan")) {
                        $X = 0.00;
                    }
                    
                    if (($ref1 eq "nan") || ($ref2 eq "nan")) {
                        $X = 0.00;
                    }
                }
                for ($n7=0;$n7<=$num_internal_nodes-1;$n7++) {
                    print (OUT2 "${'site'.$n7}\t");
                }
                print (OUT2 "$X\t");
                
                print (OUT1 "$ref_m\t$X\t");
            }
        }
        
        if ($total == 0) {
            while ($line =~ /(\w+)\s+\((\d\.\d+)\)/g) {
                $preref = $1;
                @preref = split //, $preref;
                $ref_m = @preref[$target_node-1];
                for ($n7=0;$n7<=$num_internal_nodes-1;$n7++) {
                    ${'site'.$n7} = @preref[$n7];
                }
                
                for ($n7=0;$n7<=$num_internal_nodes-1;$n7++) {
                    print (OUT2 "${'site'.$n7}\t");
                }
                print (OUT2 "0\t");
                print (OUT1 "$ref_m\t$X\t");
            }
        }
    }
    
    
    print (OUT1 "\n");
    print (OUT2 "\n");
    
}
close (OUT1);
close (OUT2);
