#merge the anc_site_probs_all_node.txt files of all species
#the output shows the ancestral inference result after the weighting for the all species and polymrophism information of each species.
$sp_num = $ARGV[0];

for ($n1=1;$n1<=$sp_num-1;$n1++) {
    
    open (IN1,"../iterative_BTWest/anc_site_probs_all_node.txt");
    @str1 = <IN1>;
    open (IN2,"../iterative_BTWest/anc_site_probs_all_node_$n1\.txt");
    @str2 = <IN2>;
    
    for ($n=0;$n<=$#str1;$n++) {
        if (@str1[$n] =~ /^(.*)\t(\w+\:.*)\n/) {
            $refA = $1;
        }
        if (@str2[$n] =~ /(\d+)\t(\d+)\t(\d+)\t(\w+\:.*)\n/) {
            $refB1 = $2;
            $refB2 = $3;
            $pattern = $4;
        }
        open (OUT1, ">>../iterative_BTWest/anc_site_probs_all_node_new.txt");
        print (OUT1 "$refA\t$refB1\t$refB2\t$pattern\n");
        close (OUT1);
    }
    unlink "../iterative_BTWest/anc_site_probs_all_node.txt";
    unlink "../iterative_BTWest/anc_site_probs_all_node_$n1\.txt";
    rename "../iterative_BTWest/anc_site_probs_all_node_new.txt", "../iterative_BTWest/anc_site_probs_all_node.txt";
}

