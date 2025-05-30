#The purpose of this code is to make the list of "current states - ancestral states configuration pairs for each site of the nucleotide sequences used in the analysis".

open (IN1,"collapse_seq.txt"); #collapse_alignment which was the input of BASEML
@str1 = <IN1>;
open (IN2,"fltrst"); #output of BASEML which shows the estimated ancestral probabilities and tree information
@str2 = <IN2>;
open (IN3,"original_aln.txt"); #original_alignment which was the input to make collapse sequences
@str3 = <IN3>;
open (IN4,"target_node"); #file showing the position of the target node (ancestor of the two collapse sequences of the focused species (Xth species in the shell script)
@str4 = <IN4>;

# read input parameters

$total_seq_num = $ARGV[0]; #number of the sequences in the collapse_alignment
$num_internal_nodes = $ARGV[1]; #number of internal nodes in the input tree of BASEML
$sample_start = $ARGV[2]; #order of the first sequence of the focused species in the original_alignemnt
$sample_end = $ARGV[3]; #order of the last sequence of the focused species in the original_alignemnt
$collapse_start = $ARGV[4]; #order of the first collapse sequence of the focused species in the collapse_alignemnt
$collapse_end = $ARGV[5]; #order of the second collapse sequence of the focused species in the collapse_alignemnt
$target_node = @str4[0];
chomp ($target_node);

#printf "$total_seq_num\t$num_internal_nodes\t$target_node\t$sample_start\t$sample_end\t$collapse_start\t$collapse_end\n";

# calculate the number of the sequences of the focused species in original_alignment
$sample_num = $sample_end - $sample_start + 1;

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
    }
    $freq_ref1 = $nref1;
    $freq_ref2 = $nref2;
    
    $n1++;
    
    #read fltrst and find the line showing the ancestral states configuration of the current states configuration of the current site
    foreach $_ (@str2) {
        #search line with the below pattern
        if ($_ =~ /(\d+)\s+$site_pattern:\s(\w+)/) {
            
            #output the info of $site, frequencies of 1st and 2nd polymorphic states (freq of the state that the 1st collapse seq has comes first)
            if (@{'SEQ'.$collapse_start}[$site] eq $ref1) {
                print (OUT2 "$site\t$freq_ref1\t$freq_ref2\t$site_pattern\:\t\t");
            }
            elsif (@{'SEQ'.$collapse_start}[$site] eq $ref2) {
                print (OUT2 "$site\t$freq_ref2\t$freq_ref1\t$site_pattern\:\t\t");
            }
            
            # read the part showing ancestral states configuration and its probability
            while ($_ =~ /(\w+)\s+\((\d\.\d+)\)/g) {
                # ancestral configuration is stored as string $preref
                $preref = $1;
                # read the probability of the configuration
                $prob = $2;
                # $preref is converted to array
                @preref = split //, $preref;
                # ref_m is the state at the ancestral node of the two collapse sequences of the focused species
                $ref_m = @preref[$target_node-1];
                # get the state of the all ancestral nodes
                for ($n7=0;$n7<=$num_internal_nodes-1;$n7++) {
                    ${'site'.$n7} = @preref[$n7];
                }
                #output the ancestral configuration and prob of the current site
                if ($preref ne "total") {
                    for ($n7=0;$n7<=$num_internal_nodes-1;$n7++) {
                        print (OUT2 "${'site'.$n7}\t");
                    }
                    print (OUT2 "$prob\t");
                }
            }
            
            print (OUT2 "\n");
            
            last;
        }
    }
}
close (OUT1);
