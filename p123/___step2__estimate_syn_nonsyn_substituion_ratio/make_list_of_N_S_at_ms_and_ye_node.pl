open (IN1,"joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_stop_codons_only_paths_with_single_or_no_change.txt");
@str1 = <IN1>;

$TC_w = 1;
$TA_w = 1;
$TG_w = 1;
$CT_w = 1;
$CA_w = 1;
$CG_w = 1;
$AT_w = 1;
$AC_w = 1;
$AG_w = 1;
$GT_w = 1;
$GC_w = 1;
$GA_w = 1;

$N_total = 0;
$S_total = 0;

foreach $_ (@str1) {
    if ($_ =~ /(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\d.*)\n/) {
        $ms = $1;
        $ye = $2;
        $m = $3;
        $s = $4;
        $prob = $5;
        
        $anc_codon = $ms;
        
        if ($anc_codon eq 'TTT') {$N = $prob*($TC_w*2+$TA_w*3+$TG_w*3);$anc_index=1;}
        if ($anc_codon eq 'TTC') {$N = $prob*($TC_w*2+$TA_w*2+$TG_w*2+$CA_w*1+$CG_w*1);$anc_index=2;}
        if ($anc_codon eq 'TTA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1);$anc_index=3;}
        if ($anc_codon eq 'TTG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*2+$GT_w*1+$GC_w*1);$anc_index=4;}
        if ($anc_codon eq 'CTT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1);$anc_index=5;}
        if ($anc_codon eq 'CTC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1);$anc_index=6;}
        if ($anc_codon eq 'CTA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CA_w*1+$CG_w*1);$anc_index=7;}
        if ($anc_codon eq 'CTG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CA_w*1+$CG_w*1);$anc_index=8;}
        if ($anc_codon eq 'ATT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*2+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=9;}
        if ($anc_codon eq 'ATC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=10;}
        if ($anc_codon eq 'ATA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*2);$anc_index=11;}
        if ($anc_codon eq 'ATG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=12;}
        if ($anc_codon eq 'GTT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=13;}
        if ($anc_codon eq 'GTC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=14;}
        if ($anc_codon eq 'GTA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=15;}
        if ($anc_codon eq 'GTG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=16;}
        if ($anc_codon eq 'TCT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1);$anc_index=17;}
        if ($anc_codon eq 'TCC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1);$anc_index=18;}
        if ($anc_codon eq 'TCA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1);$anc_index=19;}
        if ($anc_codon eq 'TCG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CG_w*1);$anc_index=20;}
        if ($anc_codon eq 'CCT') {$N = $prob*($CT_w*2+$CA_w*2+$CG_w*2);$anc_index=21;}
        if ($anc_codon eq 'CCC') {$N = $prob*($CT_w*2+$CA_w*2+$CG_w*2);$anc_index=22;}
        if ($anc_codon eq 'CCA') {$N = $prob*($CT_w*2+$CA_w*2+$CG_w*2);$anc_index=23;}
        if ($anc_codon eq 'CCG') {$N = $prob*($CT_w*2+$CA_w*2+$CG_w*2);$anc_index=24;}
        if ($anc_codon eq 'ACT') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=25;}
        if ($anc_codon eq 'ACC') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=26;}
        if ($anc_codon eq 'ACA') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=27;}
        if ($anc_codon eq 'ACG') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=28;}
        if ($anc_codon eq 'GCT') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=29;}
        if ($anc_codon eq 'GCC') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=30;}
        if ($anc_codon eq 'GCA') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=31;}
        if ($anc_codon eq 'GCG') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=32;}
        if ($anc_codon eq 'TAT') {$N = $prob*($TC_w*1+$TA_w*2+$TG_w*2+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=33;}
        if ($anc_codon eq 'TAC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=34;}
        if ($anc_codon eq 'TAA') {$N = $prob*0;$anc_index=35;}
        if ($anc_codon eq 'TAG') {$N = $prob*0;$anc_index=36;}
        if ($anc_codon eq 'CAT') {$N = $prob*($TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=37;}
        if ($anc_codon eq 'CAC') {$N = $prob*($CT_w*1+$CA_w*2+$CG_w*2+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=38;}
        if ($anc_codon eq 'CAA') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*2+$AC_w*2+$AG_w*1);$anc_index=39;}
        if ($anc_codon eq 'CAG') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1);$anc_index=40;}
        if ($anc_codon eq 'AAT') {$N = $prob*($TA_w*1+$TG_w*1+$AT_w*2+$AC_w*2+$AG_w*2);$anc_index=41;}
        if ($anc_codon eq 'AAC') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*2+$AC_w*2+$AG_w*2);$anc_index=42;}
        if ($anc_codon eq 'AAA') {$N = $prob*($AT_w*2+$AC_w*3+$AG_w*2);$anc_index=43;}
        if ($anc_codon eq 'AAG') {$N = $prob*($AT_w*1+$AC_w*2+$AG_w*2+$GT_w*1+$GC_w*1);$anc_index=44;}
        if ($anc_codon eq 'GAT') {$N = $prob*($TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=45;}
        if ($anc_codon eq 'GAC') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=46;}
        if ($anc_codon eq 'GAA') {$N = $prob*($AT_w*2+$AC_w*2+$AG_w*1+$GC_w*1+$GA_w*1);$anc_index=47;}
        if ($anc_codon eq 'GAG') {$N = $prob*($AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*2+$GA_w*1);$anc_index=48;}
        if ($anc_codon eq 'TGT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*2+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=49;}
        if ($anc_codon eq 'TGC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=50;}
        if ($anc_codon eq 'TGA') {$N = $prob*0;$anc_index=51;}
        if ($anc_codon eq 'TGG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*2+$GC_w*2);$anc_index=52;}
        if ($anc_codon eq 'CGT') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=53;}
        if ($anc_codon eq 'CGC') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=54;}
        if ($anc_codon eq 'CGA') {$N = $prob*($CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=55;}
        if ($anc_codon eq 'CGG') {$N = $prob*($CT_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=56;}
        if ($anc_codon eq 'AGT') {$N = $prob*($TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=57;}
        if ($anc_codon eq 'AGC') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=58;}
        if ($anc_codon eq 'AGA') {$N = $prob*($AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=59;}
        if ($anc_codon eq 'AGG') {$N = $prob*($AT_w*1+$AG_w*1+$GT_w*2+$GC_w*2+$GA_w*1);$anc_index=60;}
        if ($anc_codon eq 'GGT') {$N = $prob*($GT_w*2+$GC_w*2+$GA_w*2);$anc_index=61;}
        if ($anc_codon eq 'GGC') {$N = $prob*($GT_w*2+$GC_w*2+$GA_w*2);$anc_index=62;}
        if ($anc_codon eq 'GGA') {$N = $prob*($GT_w*1+$GC_w*2+$GA_w*2);$anc_index=63;}
        if ($anc_codon eq 'GGG') {$N = $prob*($GT_w*2+$GC_w*2+$GA_w*2);$anc_index=64;}
        
        if ($anc_codon eq 'TTT') {$S = $prob*1;$anc_index=1;}
        if ($anc_codon eq 'TTC') {$S = $prob*1;$anc_index=2;}
        if ($anc_codon eq 'TTA') {$S = $prob*2;$anc_index=3;}
        if ($anc_codon eq 'TTG') {$S = $prob*2;$anc_index=4;}
        if ($anc_codon eq 'CTT') {$S = $prob*3;$anc_index=5;}
        if ($anc_codon eq 'CTC') {$S = $prob*3;$anc_index=6;}
        if ($anc_codon eq 'CTA') {$S = $prob*4;$anc_index=7;}
        if ($anc_codon eq 'CTG') {$S = $prob*4;$anc_index=8;}
        if ($anc_codon eq 'ATT') {$S = $prob*2;$anc_index=9;}
        if ($anc_codon eq 'ATC') {$S = $prob*2;$anc_index=10;}
        if ($anc_codon eq 'ATA') {$S = $prob*2;$anc_index=11;}
        if ($anc_codon eq 'ATG') {$S = $prob*0;$anc_index=12;}
        if ($anc_codon eq 'GTT') {$S = $prob*3;$anc_index=13;}
        if ($anc_codon eq 'GTC') {$S = $prob*3;$anc_index=14;}
        if ($anc_codon eq 'GTA') {$S = $prob*3;$anc_index=15;}
        if ($anc_codon eq 'GTG') {$S = $prob*3;$anc_index=16;}
        if ($anc_codon eq 'TCT') {$S = $prob*3;$anc_index=17;}
        if ($anc_codon eq 'TCC') {$S = $prob*3;$anc_index=18;}
        if ($anc_codon eq 'TCA') {$S = $prob*3;$anc_index=19;}
        if ($anc_codon eq 'TCG') {$S = $prob*3;$anc_index=20;}
        if ($anc_codon eq 'CCT') {$S = $prob*3;$anc_index=21;}
        if ($anc_codon eq 'CCC') {$S = $prob*3;$anc_index=22;}
        if ($anc_codon eq 'CCA') {$S = $prob*3;$anc_index=23;}
        if ($anc_codon eq 'CCG') {$S = $prob*3;$anc_index=24;}
        if ($anc_codon eq 'ACT') {$S = $prob*3;$anc_index=25;}
        if ($anc_codon eq 'ACC') {$S = $prob*3;$anc_index=26;}
        if ($anc_codon eq 'ACA') {$S = $prob*3;$anc_index=27;}
        if ($anc_codon eq 'ACG') {$S = $prob*3;$anc_index=28;}
        if ($anc_codon eq 'GCT') {$S = $prob*3;$anc_index=29;}
        if ($anc_codon eq 'GCC') {$S = $prob*3;$anc_index=30;}
        if ($anc_codon eq 'GCA') {$S = $prob*3;$anc_index=31;}
        if ($anc_codon eq 'GCG') {$S = $prob*3;$anc_index=32;}
        if ($anc_codon eq 'TAT') {$S = $prob*1;$anc_index=33;}
        if ($anc_codon eq 'TAC') {$S = $prob*1;$anc_index=34;}
        if ($anc_codon eq 'TAA') {$S = $prob*0;$anc_index=35;}
        if ($anc_codon eq 'TAG') {$S = $prob*0;$anc_index=36;}
        if ($anc_codon eq 'CAT') {$S = $prob*1;$anc_index=37;}
        if ($anc_codon eq 'CAC') {$S = $prob*1;$anc_index=38;}
        if ($anc_codon eq 'CAA') {$S = $prob*1;$anc_index=39;}
        if ($anc_codon eq 'CAG') {$S = $prob*1;$anc_index=40;}
        if ($anc_codon eq 'AAT') {$S = $prob*1;$anc_index=41;}
        if ($anc_codon eq 'AAC') {$S = $prob*1;$anc_index=42;}
        if ($anc_codon eq 'AAA') {$S = $prob*1;$anc_index=43;}
        if ($anc_codon eq 'AAG') {$S = $prob*1;$anc_index=44;}
        if ($anc_codon eq 'GAT') {$S = $prob*1;$anc_index=45;}
        if ($anc_codon eq 'GAC') {$S = $prob*1;$anc_index=46;}
        if ($anc_codon eq 'GAA') {$S = $prob*1;$anc_index=47;}
        if ($anc_codon eq 'GAG') {$S = $prob*1;$anc_index=48;}
        if ($anc_codon eq 'TGT') {$S = $prob*1;$anc_index=49;}
        if ($anc_codon eq 'TGC') {$S = $prob*1;$anc_index=50;}
        if ($anc_codon eq 'TGA') {$S = $prob*0;$anc_index=51;}
        if ($anc_codon eq 'TGG') {$S = $prob*0;$anc_index=52;}
        if ($anc_codon eq 'CGT') {$S = $prob*3;$anc_index=53;}
        if ($anc_codon eq 'CGC') {$S = $prob*3;$anc_index=54;}
        if ($anc_codon eq 'CGA') {$S = $prob*4;$anc_index=55;}
        if ($anc_codon eq 'CGG') {$S = $prob*4;$anc_index=56;}
        if ($anc_codon eq 'AGT') {$S = $prob*1;$anc_index=57;}
        if ($anc_codon eq 'AGC') {$S = $prob*1;$anc_index=58;}
        if ($anc_codon eq 'AGA') {$S = $prob*2;$anc_index=59;}
        if ($anc_codon eq 'AGG') {$S = $prob*2;$anc_index=60;}
        if ($anc_codon eq 'GGT') {$S = $prob*3;$anc_index=61;}
        if ($anc_codon eq 'GGC') {$S = $prob*3;$anc_index=62;}
        if ($anc_codon eq 'GGA') {$S = $prob*3;$anc_index=63;}
        if ($anc_codon eq 'GGG') {$S = $prob*3;$anc_index=64;}
        
        $N_total = $N_total + $N;
        $S_total = $S_total + $S;
    }
}
printf ("ms node\tN=$N_total\tS=$S_total\n");


$N_total = 0;
$S_total = 0;

foreach $_ (@str1) {
    if ($_ =~ /(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\w\w\w)\t(\d.*)\n/) {
        $ms = $1;
        $ye = $2;
        $m = $3;
        $s = $4;
        $prob = $5;
        
        $anc_codon = $ye;
        
        if ($anc_codon eq 'TTT') {$N = $prob*($TC_w*2+$TA_w*3+$TG_w*3);$anc_index=1;}
        if ($anc_codon eq 'TTC') {$N = $prob*($TC_w*2+$TA_w*2+$TG_w*2+$CA_w*1+$CG_w*1);$anc_index=2;}
        if ($anc_codon eq 'TTA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1);$anc_index=3;}
        if ($anc_codon eq 'TTG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*2+$GT_w*1+$GC_w*1);$anc_index=4;}
        if ($anc_codon eq 'CTT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1);$anc_index=5;}
        if ($anc_codon eq 'CTC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1);$anc_index=6;}
        if ($anc_codon eq 'CTA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CA_w*1+$CG_w*1);$anc_index=7;}
        if ($anc_codon eq 'CTG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CA_w*1+$CG_w*1);$anc_index=8;}
        if ($anc_codon eq 'ATT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*2+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=9;}
        if ($anc_codon eq 'ATC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=10;}
        if ($anc_codon eq 'ATA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*2);$anc_index=11;}
        if ($anc_codon eq 'ATG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=12;}
        if ($anc_codon eq 'GTT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=13;}
        if ($anc_codon eq 'GTC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=14;}
        if ($anc_codon eq 'GTA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=15;}
        if ($anc_codon eq 'GTG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=16;}
        if ($anc_codon eq 'TCT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1);$anc_index=17;}
        if ($anc_codon eq 'TCC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1);$anc_index=18;}
        if ($anc_codon eq 'TCA') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1);$anc_index=19;}
        if ($anc_codon eq 'TCG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CT_w*1+$CG_w*1);$anc_index=20;}
        if ($anc_codon eq 'CCT') {$N = $prob*($CT_w*2+$CA_w*2+$CG_w*2);$anc_index=21;}
        if ($anc_codon eq 'CCC') {$N = $prob*($CT_w*2+$CA_w*2+$CG_w*2);$anc_index=22;}
        if ($anc_codon eq 'CCA') {$N = $prob*($CT_w*2+$CA_w*2+$CG_w*2);$anc_index=23;}
        if ($anc_codon eq 'CCG') {$N = $prob*($CT_w*2+$CA_w*2+$CG_w*2);$anc_index=24;}
        if ($anc_codon eq 'ACT') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=25;}
        if ($anc_codon eq 'ACC') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=26;}
        if ($anc_codon eq 'ACA') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=27;}
        if ($anc_codon eq 'ACG') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=28;}
        if ($anc_codon eq 'GCT') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=29;}
        if ($anc_codon eq 'GCC') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=30;}
        if ($anc_codon eq 'GCA') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=31;}
        if ($anc_codon eq 'GCG') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=32;}
        if ($anc_codon eq 'TAT') {$N = $prob*($TC_w*1+$TA_w*2+$TG_w*2+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=33;}
        if ($anc_codon eq 'TAC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=34;}
        if ($anc_codon eq 'TAA') {$N = $prob*0;$anc_index=35;}
        if ($anc_codon eq 'TAG') {$N = $prob*0;$anc_index=36;}
        if ($anc_codon eq 'CAT') {$N = $prob*($TA_w*1+$TG_w*1+$CT_w*1+$CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=37;}
        if ($anc_codon eq 'CAC') {$N = $prob*($CT_w*1+$CA_w*2+$CG_w*2+$AT_w*1+$AC_w*1+$AG_w*1);$anc_index=38;}
        if ($anc_codon eq 'CAA') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*2+$AC_w*2+$AG_w*1);$anc_index=39;}
        if ($anc_codon eq 'CAG') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1);$anc_index=40;}
        if ($anc_codon eq 'AAT') {$N = $prob*($TA_w*1+$TG_w*1+$AT_w*2+$AC_w*2+$AG_w*2);$anc_index=41;}
        if ($anc_codon eq 'AAC') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*2+$AC_w*2+$AG_w*2);$anc_index=42;}
        if ($anc_codon eq 'AAA') {$N = $prob*($AT_w*2+$AC_w*3+$AG_w*2);$anc_index=43;}
        if ($anc_codon eq 'AAG') {$N = $prob*($AT_w*1+$AC_w*2+$AG_w*2+$GT_w*1+$GC_w*1);$anc_index=44;}
        if ($anc_codon eq 'GAT') {$N = $prob*($TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=45;}
        if ($anc_codon eq 'GAC') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=46;}
        if ($anc_codon eq 'GAA') {$N = $prob*($AT_w*2+$AC_w*2+$AG_w*1+$GC_w*1+$GA_w*1);$anc_index=47;}
        if ($anc_codon eq 'GAG') {$N = $prob*($AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*2+$GA_w*1);$anc_index=48;}
        if ($anc_codon eq 'TGT') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*2+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=49;}
        if ($anc_codon eq 'TGC') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=50;}
        if ($anc_codon eq 'TGA') {$N = $prob*0;$anc_index=51;}
        if ($anc_codon eq 'TGG') {$N = $prob*($TC_w*1+$TA_w*1+$TG_w*1+$GT_w*2+$GC_w*2);$anc_index=52;}
        if ($anc_codon eq 'CGT') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=53;}
        if ($anc_codon eq 'CGC') {$N = $prob*($CT_w*1+$CA_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=54;}
        if ($anc_codon eq 'CGA') {$N = $prob*($CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=55;}
        if ($anc_codon eq 'CGG') {$N = $prob*($CT_w*1+$CG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=56;}
        if ($anc_codon eq 'AGT') {$N = $prob*($TA_w*1+$TG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=57;}
        if ($anc_codon eq 'AGC') {$N = $prob*($CA_w*1+$CG_w*1+$AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=58;}
        if ($anc_codon eq 'AGA') {$N = $prob*($AT_w*1+$AC_w*1+$AG_w*1+$GT_w*1+$GC_w*1+$GA_w*1);$anc_index=59;}
        if ($anc_codon eq 'AGG') {$N = $prob*($AT_w*1+$AG_w*1+$GT_w*2+$GC_w*2+$GA_w*1);$anc_index=60;}
        if ($anc_codon eq 'GGT') {$N = $prob*($GT_w*2+$GC_w*2+$GA_w*2);$anc_index=61;}
        if ($anc_codon eq 'GGC') {$N = $prob*($GT_w*2+$GC_w*2+$GA_w*2);$anc_index=62;}
        if ($anc_codon eq 'GGA') {$N = $prob*($GT_w*1+$GC_w*2+$GA_w*2);$anc_index=63;}
        if ($anc_codon eq 'GGG') {$N = $prob*($GT_w*2+$GC_w*2+$GA_w*2);$anc_index=64;}
        
        if ($anc_codon eq 'TTT') {$S = $prob*1;$anc_index=1;}
        if ($anc_codon eq 'TTC') {$S = $prob*1;$anc_index=2;}
        if ($anc_codon eq 'TTA') {$S = $prob*2;$anc_index=3;}
        if ($anc_codon eq 'TTG') {$S = $prob*2;$anc_index=4;}
        if ($anc_codon eq 'CTT') {$S = $prob*3;$anc_index=5;}
        if ($anc_codon eq 'CTC') {$S = $prob*3;$anc_index=6;}
        if ($anc_codon eq 'CTA') {$S = $prob*4;$anc_index=7;}
        if ($anc_codon eq 'CTG') {$S = $prob*4;$anc_index=8;}
        if ($anc_codon eq 'ATT') {$S = $prob*2;$anc_index=9;}
        if ($anc_codon eq 'ATC') {$S = $prob*2;$anc_index=10;}
        if ($anc_codon eq 'ATA') {$S = $prob*2;$anc_index=11;}
        if ($anc_codon eq 'ATG') {$S = $prob*0;$anc_index=12;}
        if ($anc_codon eq 'GTT') {$S = $prob*3;$anc_index=13;}
        if ($anc_codon eq 'GTC') {$S = $prob*3;$anc_index=14;}
        if ($anc_codon eq 'GTA') {$S = $prob*3;$anc_index=15;}
        if ($anc_codon eq 'GTG') {$S = $prob*3;$anc_index=16;}
        if ($anc_codon eq 'TCT') {$S = $prob*3;$anc_index=17;}
        if ($anc_codon eq 'TCC') {$S = $prob*3;$anc_index=18;}
        if ($anc_codon eq 'TCA') {$S = $prob*3;$anc_index=19;}
        if ($anc_codon eq 'TCG') {$S = $prob*3;$anc_index=20;}
        if ($anc_codon eq 'CCT') {$S = $prob*3;$anc_index=21;}
        if ($anc_codon eq 'CCC') {$S = $prob*3;$anc_index=22;}
        if ($anc_codon eq 'CCA') {$S = $prob*3;$anc_index=23;}
        if ($anc_codon eq 'CCG') {$S = $prob*3;$anc_index=24;}
        if ($anc_codon eq 'ACT') {$S = $prob*3;$anc_index=25;}
        if ($anc_codon eq 'ACC') {$S = $prob*3;$anc_index=26;}
        if ($anc_codon eq 'ACA') {$S = $prob*3;$anc_index=27;}
        if ($anc_codon eq 'ACG') {$S = $prob*3;$anc_index=28;}
        if ($anc_codon eq 'GCT') {$S = $prob*3;$anc_index=29;}
        if ($anc_codon eq 'GCC') {$S = $prob*3;$anc_index=30;}
        if ($anc_codon eq 'GCA') {$S = $prob*3;$anc_index=31;}
        if ($anc_codon eq 'GCG') {$S = $prob*3;$anc_index=32;}
        if ($anc_codon eq 'TAT') {$S = $prob*1;$anc_index=33;}
        if ($anc_codon eq 'TAC') {$S = $prob*1;$anc_index=34;}
        if ($anc_codon eq 'TAA') {$S = $prob*0;$anc_index=35;}
        if ($anc_codon eq 'TAG') {$S = $prob*0;$anc_index=36;}
        if ($anc_codon eq 'CAT') {$S = $prob*1;$anc_index=37;}
        if ($anc_codon eq 'CAC') {$S = $prob*1;$anc_index=38;}
        if ($anc_codon eq 'CAA') {$S = $prob*1;$anc_index=39;}
        if ($anc_codon eq 'CAG') {$S = $prob*1;$anc_index=40;}
        if ($anc_codon eq 'AAT') {$S = $prob*1;$anc_index=41;}
        if ($anc_codon eq 'AAC') {$S = $prob*1;$anc_index=42;}
        if ($anc_codon eq 'AAA') {$S = $prob*1;$anc_index=43;}
        if ($anc_codon eq 'AAG') {$S = $prob*1;$anc_index=44;}
        if ($anc_codon eq 'GAT') {$S = $prob*1;$anc_index=45;}
        if ($anc_codon eq 'GAC') {$S = $prob*1;$anc_index=46;}
        if ($anc_codon eq 'GAA') {$S = $prob*1;$anc_index=47;}
        if ($anc_codon eq 'GAG') {$S = $prob*1;$anc_index=48;}
        if ($anc_codon eq 'TGT') {$S = $prob*1;$anc_index=49;}
        if ($anc_codon eq 'TGC') {$S = $prob*1;$anc_index=50;}
        if ($anc_codon eq 'TGA') {$S = $prob*0;$anc_index=51;}
        if ($anc_codon eq 'TGG') {$S = $prob*0;$anc_index=52;}
        if ($anc_codon eq 'CGT') {$S = $prob*3;$anc_index=53;}
        if ($anc_codon eq 'CGC') {$S = $prob*3;$anc_index=54;}
        if ($anc_codon eq 'CGA') {$S = $prob*4;$anc_index=55;}
        if ($anc_codon eq 'CGG') {$S = $prob*4;$anc_index=56;}
        if ($anc_codon eq 'AGT') {$S = $prob*1;$anc_index=57;}
        if ($anc_codon eq 'AGC') {$S = $prob*1;$anc_index=58;}
        if ($anc_codon eq 'AGA') {$S = $prob*2;$anc_index=59;}
        if ($anc_codon eq 'AGG') {$S = $prob*2;$anc_index=60;}
        if ($anc_codon eq 'GGT') {$S = $prob*3;$anc_index=61;}
        if ($anc_codon eq 'GGC') {$S = $prob*3;$anc_index=62;}
        if ($anc_codon eq 'GGA') {$S = $prob*3;$anc_index=63;}
        if ($anc_codon eq 'GGG') {$S = $prob*3;$anc_index=64;}
        
        $N_total = $N_total + $N;
        $S_total = $S_total + $S;
    }
}

printf ("ye node\tN=$N_total\tS=$S_total\n");

