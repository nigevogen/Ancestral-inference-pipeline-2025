#make fltrst like file which keep the fltrst format but list the information of each site in one line

open (IN1,"anc_site_probs_all_node.txt");
@str1 = <IN1>;

$site = 0;
foreach $_ (@str1) {
    if ($_ =~ /(\d+)\t(\d+)\t(\d+)\t(\w+)\:\t\t(.*)\n/) {
        
        $ref1 = $2;
        $ref2 = $3;
        $pattern = $4;
        $line = $5;
        @line = split /\t/, $line;
        
        $count = 1;
        
        open (OUT1, ">>fltrst_re");
        print (OUT1 "$site     $count     $pattern\: ");
        
        foreach $_ (@line) {
            if ($_ =~ /^(\D)$/) {
                print (OUT1 "$1");
            }
            if ($_ =~ /^(\d\.\d+)$/) {
                print (OUT1 " \($1\) ");
            }
            if ($_ =~ /^(1)$/) {
                print (OUT1 " \(1.000\) ");
            }
            if ($_ =~ /^(0)$/) {
                print (OUT1 " \(0.000\) ");
            }
        }
        print (OUT1 "\n");
        close (OUT1);
        $site++;
    }
}
