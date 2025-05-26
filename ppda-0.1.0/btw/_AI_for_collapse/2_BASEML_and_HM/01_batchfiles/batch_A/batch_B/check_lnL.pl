#make INF_PATTERNS file for each replicate
open (IN1,"mlf");
open (IN2,"highestL");
my @str1 = <IN1>;
my @str2 = <IN2>;

$hlnL=-9999999999;
foreach $_ (@str2) {
    if ($_ =~ /lnL=(\-\d+\.\d+)/) {
        $hlnL = $1;
    }
}

foreach $_ (@str1) {
    if ($_ =~ /lnL\(ntime:.*\):\s+(\-\d+\.\d+)/) {
        if ($1 > $hlnL) {
            print "1\n";
            open (OUT1, ">highestL");
            print (OUT1 "lnL=$1\n");
            close (OUT1);
        }
        else {
            print "0\n";
        }
    }
}

