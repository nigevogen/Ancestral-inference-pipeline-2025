open (IN1,"fltrst");
@str1 = <IN1>;

$num_seq = $ARGV[0]; #number of sequences
$pos_collapse_A = $ARGV[1]; #order of the 1st collapse sequence in the collapse_alingment
$pos_collapse_B = $ARGV[2]; #order of the 2nd collapse sequence in the collapse_alingment

#search the line showing the relationship of ancestral-derived node in flrst
#and find the common ancestral node of the two collapse sequence
foreach $_ (@str1) {
    if ($_ =~ /\s\s(\d+)\.\.($pos_collapse_A)\s/) {
        $anc_node_A = $1;
    }
    if ($_ =~ /\s\s(\d+)\.\.($pos_collapse_B)\s/) {
        $anc_node_B = $1;
    }
}

#if there is a common ancestral node, store its order among the all ancestral nodes in the output file
if ($anc_node_A == $anc_node_B) {
    #each node is named as 1, 2, 3, ... from derived to ancestral nodes. So the order of the common ancestral node among the all ancestral nodes is calculated as below
    $target_node = $anc_node_A - $num_seq;
    open (OUT1, ">>target_node");
    print (OUT1 "$target_node\n");
    close (OUT1);
}

#if there is no common ancestral node, output the error message
if ($anc_node_A != $anc_node_B) {
    printf "ERROR: could not find the ancestral node of the collapse sequences in fltrst\n";
}
