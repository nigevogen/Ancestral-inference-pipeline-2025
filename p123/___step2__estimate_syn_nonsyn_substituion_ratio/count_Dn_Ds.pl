open (IN1,"list_of_codon_change_with_probability_S_N_Dere_lineage_for_NG_dNdS_calculation.txt");
@str1 = <IN1>;

$S = 0;
$N = 0;

foreach $_ (@str1) {
    @line = split /\t/, $_;
    
    if (@line[6] eq 'S') {
        $S = $S + @line[5];
    }
    if (@line[6] eq 'N') {
        $N = $N + @line[5];
    }
}

print "Dere lineage Dn=$N\tDs=$S\n";

open (IN1,"list_of_codon_change_with_probability_S_N_Dyak_lineage_for_NG_dNdS_calculation.txt");
@str1 = <IN1>;

$S = 0;
$N = 0;

foreach $_ (@str1) {
    @line = split /\t/, $_;
    
    if (@line[6] eq 'S') {
        $S = $S + @line[5];
    }
    if (@line[6] eq 'N') {
        $N = $N + @line[5];
    }
}

print "Dyak lineage Dn=$N\tDs=$S\n";

open (IN1,"list_of_codon_change_with_probability_S_N_Dmel_lineage_for_NG_dNdS_calculation.txt");
@str1 = <IN1>;

$S = 0;
$N = 0;

foreach $_ (@str1) {
    @line = split /\t/, $_;
    
    if (@line[6] eq 'S') {
        $S = $S + @line[5];
    }
    if (@line[6] eq 'N') {
        $N = $N + @line[5];
    }
}

print "Dmel lineage Dn=$N\tDs=$S\n";

open (IN1,"list_of_codon_change_with_probability_S_N_Dsim_lineage_for_NG_dNdS_calculation.txt");
@str1 = <IN1>;

$S = 0;
$N = 0;

foreach $_ (@str1) {
    @line = split /\t/, $_;
    
    if (@line[6] eq 'S') {
        $S = $S + @line[5];
    }
    if (@line[6] eq 'N') {
        $N = $N + @line[5];
    }
}

print "Dsim lineage Dn=$N\tDs=$S\n";

open (IN1,"list_of_codon_change_with_probability_S_N_ye-ms_lineage_for_NG_dNdS_calculation.txt");
@str1 = <IN1>;

$S = 0;
$N = 0;

foreach $_ (@str1) {
    @line = split /\t/, $_;
    
    if (@line[6] eq 'S') {
        $S = $S + @line[5];
    }
    if (@line[6] eq 'N') {
        $N = $N + @line[5];
    }
}

print "ye-ms lineage Dn=$N\tDs=$S\n";
