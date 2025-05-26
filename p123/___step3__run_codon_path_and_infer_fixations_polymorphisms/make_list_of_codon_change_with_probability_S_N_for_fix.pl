open (IN1,"list_of_codon_change_with_probability_Dmel_lineage.txt");
@str1 = <IN1>;

foreach $_ (@str1) {
    if ($_ =~ /(\d+)\t(\w\w\w)\t(\w\w\w)\t(\d.*)\n/) {
        $codon_pos = $1;
        $anc_codon = $2;
        $der_codon = $3;
        $prob = $4;
        
        $type1 = '';
        $type2 = '';
        
        if (($anc_codon eq 'TTT')||($anc_codon eq 'TTC')) {$anc_AA = 'PHE';}
        if (($anc_codon eq 'TTA')||($anc_codon eq 'TTG')) {$anc_AA = 'LEUa';}
        if (($anc_codon eq 'CTT')||($anc_codon eq 'CTC')||($anc_codon eq 'CTA')||($anc_codon eq 'CTG')) {$anc_AA = 'LEUb';}
        if (($anc_codon eq 'ATT')||($anc_codon eq 'ATC')||($anc_codon eq 'ATA')) {$anc_AA = 'ILE';}
        if (($anc_codon eq 'ATG')) {$anc_AA = 'MET';}
        if (($anc_codon eq 'GTT')||($anc_codon eq 'GTC')||($anc_codon eq 'GTA')||($anc_codon eq 'GTG')) {$anc_AA = 'VAL';}
        if (($anc_codon eq 'TCT')||($anc_codon eq 'TCC')||($anc_codon eq 'TCA')||($anc_codon eq 'TCG')) {$anc_AA = 'SERb';}
        if (($anc_codon eq 'CCT')||($anc_codon eq 'CCC')||($anc_codon eq 'CCA')||($anc_codon eq 'CCG')) {$anc_AA = 'PRO';}
        if (($anc_codon eq 'ACT')||($anc_codon eq 'ACC')||($anc_codon eq 'ACA')||($anc_codon eq 'ACG')) {$anc_AA = 'THR';}
        if (($anc_codon eq 'GCT')||($anc_codon eq 'GCC')||($anc_codon eq 'GCA')||($anc_codon eq 'GCG')) {$anc_AA = 'ALA';}
        if (($anc_codon eq 'TAT')||($anc_codon eq 'TAC')) {$anc_AA = 'TYR';}
        if (($anc_codon eq 'TAA')||($anc_codon eq 'TAG')) {$anc_AA = 'stop';}
        if (($anc_codon eq 'CAT')||($anc_codon eq 'CAC')) {$anc_AA = 'HIS';}
        if (($anc_codon eq 'CAA')||($anc_codon eq 'CAG')) {$anc_AA = 'GLN';}
        if (($anc_codon eq 'AAT')||($anc_codon eq 'AAC')) {$anc_AA = 'ASN';}
        if (($anc_codon eq 'AAA')||($anc_codon eq 'AAG')) {$anc_AA = 'LYS';}
        if (($anc_codon eq 'GAT')||($anc_codon eq 'GAC')) {$anc_AA = 'ASP';}
        if (($anc_codon eq 'GAA')||($anc_codon eq 'GAG')) {$anc_AA = 'GLU';}
        if (($anc_codon eq 'TGT')||($anc_codon eq 'TGC')) {$anc_AA = 'CYS';}
        if (($anc_codon eq 'TGA')) {$anc_AA = 'stop';}
        if (($anc_codon eq 'TGG')) {$anc_AA = 'TRP';}
        if (($anc_codon eq 'CGT')||($anc_codon eq 'CGC')||($anc_codon eq 'CGA')||($anc_codon eq 'CGG')) {$anc_AA = 'ARGb';}
        if (($anc_codon eq 'AGT')||($anc_codon eq 'AGC')) {$anc_AA = 'SERa';}
        if (($anc_codon eq 'AGA')||($anc_codon eq 'AGG')) {$anc_AA = 'ARGa';}
        if (($anc_codon eq 'GGT')||($anc_codon eq 'GGC')||($anc_codon eq 'GGA')||($anc_codon eq 'GGG')) {$anc_AA = 'GLY';}
        
        
        if (($der_codon eq 'TTT')||($der_codon eq 'TTC')) {$der_AA = 'PHE';}
        if (($der_codon eq 'TTA')||($der_codon eq 'TTG')) {$der_AA = 'LEUa';}
        if (($der_codon eq 'CTT')||($der_codon eq 'CTC')||($der_codon eq 'CTA')||($der_codon eq 'CTG')) {$der_AA = 'LEUb';}
        if (($der_codon eq 'ATT')||($der_codon eq 'ATC')||($der_codon eq 'ATA')) {$der_AA = 'ILE';}
        if (($der_codon eq 'ATG')) {$der_AA = 'MET';}
        if (($der_codon eq 'GTT')||($der_codon eq 'GTC')||($der_codon eq 'GTA')||($der_codon eq 'GTG')) {$der_AA = 'VAL';}
        if (($der_codon eq 'TCT')||($der_codon eq 'TCC')||($der_codon eq 'TCA')||($der_codon eq 'TCG')) {$der_AA = 'SERb';}
        if (($der_codon eq 'CCT')||($der_codon eq 'CCC')||($der_codon eq 'CCA')||($der_codon eq 'CCG')) {$der_AA = 'PRO';}
        if (($der_codon eq 'ACT')||($der_codon eq 'ACC')||($der_codon eq 'ACA')||($der_codon eq 'ACG')) {$der_AA = 'THR';}
        if (($der_codon eq 'GCT')||($der_codon eq 'GCC')||($der_codon eq 'GCA')||($der_codon eq 'GCG')) {$der_AA = 'ALA';}
        if (($der_codon eq 'TAT')||($der_codon eq 'TAC')) {$der_AA = 'TYR';}
        if (($der_codon eq 'TAA')||($der_codon eq 'TAG')) {$der_AA = 'stop';}
        if (($der_codon eq 'CAT')||($der_codon eq 'CAC')) {$der_AA = 'HIS';}
        if (($der_codon eq 'CAA')||($der_codon eq 'CAG')) {$der_AA = 'GLN';}
        if (($der_codon eq 'AAT')||($der_codon eq 'AAC')) {$der_AA = 'ASN';}
        if (($der_codon eq 'AAA')||($der_codon eq 'AAG')) {$der_AA = 'LYS';}
        if (($der_codon eq 'GAT')||($der_codon eq 'GAC')) {$der_AA = 'ASP';}
        if (($der_codon eq 'GAA')||($der_codon eq 'GAG')) {$der_AA = 'GLU';}
        if (($der_codon eq 'TGT')||($der_codon eq 'TGC')) {$der_AA = 'CYS';}
        if (($der_codon eq 'TGA')) {$der_AA = 'stop';}
        if (($der_codon eq 'TGG')) {$der_AA = 'TRP';}
        if (($der_codon eq 'CGT')||($der_codon eq 'CGC')||($der_codon eq 'CGA')||($der_codon eq 'CGG')) {$der_AA = 'ARGb';}
        if (($der_codon eq 'AGT')||($der_codon eq 'AGC')) {$der_AA = 'SERa';}
        if (($der_codon eq 'AGA')||($der_codon eq 'AGG')) {$der_AA = 'ARGa';}
        if (($der_codon eq 'GGT')||($der_codon eq 'GGC')||($der_codon eq 'GGA')||($der_codon eq 'GGG')) {$der_AA = 'GLY';}
        
        if ($anc_AA ne $der_AA) {
            if (($anc_AA =~ /LEU/)&&($der_AA =~ /LEU/)) {
                $type1 = 'S';
                $type2 = 'non 3rd pos change';
            }
            elsif (($anc_AA =~ /SER/)&&($der_AA =~ /SER/)) {
                $type1 = 'S';
                $type2 = 'non 3rd pos change';
            }
            elsif (($anc_AA =~ /ARG/)&&($der_AA =~ /ARG/)) {
                $type1 = 'S';
                $type2 = 'non 3rd pos change';
            }
            else {
                $type1 = 'N';
                $type2 = '';
            }
        }
        if ($anc_AA eq $der_AA) {
            $type1 = 'S';
            if ($anc_codon =~ /\w\wT/) {
                if ($der_codon =~ /\w\wC/) {
                    $type2 = 'T->C at 3rd pos';
                }
                if ($der_codon =~ /\w\wA/) {
                    $type2 = 'T->A at 3rd pos';
                }
                if ($der_codon =~ /\w\wG/) {
                    $type2 = 'T->G at 3rd pos';
                }
            }
            if ($anc_codon =~ /\w\wC/) {
                if ($der_codon =~ /\w\wT/) {
                    $type2 = 'C->T at 3rd pos';
                }
                if ($der_codon =~ /\w\wA/) {
                    $type2 = 'C->A at 3rd pos';
                }
                if ($der_codon =~ /\w\wG/) {
                    $type2 = 'C->G at 3rd pos';
                }
            }
            if ($anc_codon =~ /\w\wA/) {
                if ($der_codon =~ /\w\wT/) {
                    $type2 = 'A->T at 3rd pos';
                }
                if ($der_codon =~ /\w\wC/) {
                    $type2 = 'A->C at 3rd pos';
                }
                if ($der_codon =~ /\w\wG/) {
                    $type2 = 'A->G at 3rd pos';
                }
            }
            if ($anc_codon =~ /\w\wG/) {
                if ($der_codon =~ /\w\wT/) {
                    $type2 = 'G->T at 3rd pos';
                }
                if ($der_codon =~ /\w\wC/) {
                    $type2 = 'G->C at 3rd pos';
                }
                if ($der_codon =~ /\w\wA/) {
                    $type2 = 'G->A at 3rd pos';
                }
            }
        }
        
        open (OUT1, ">>list_of_codon_change_with_probability_S_N_Dmel_lineage.txt");
        print (OUT1 "$codon_pos\t$anc_codon\t$der_codon\t$anc_AA\t$der_AA\t$prob\t$type1\t$type2\n");
        close (OUT1);
    }
}

open (IN1,"list_of_codon_change_with_probability_Dsim_lineage.txt");
@str1 = <IN1>;

foreach $_ (@str1) {
    if ($_ =~ /(\d+)\t(\w\w\w)\t(\w\w\w)\t(\d.*)\n/) {
        $codon_pos = $1;
        $anc_codon = $2;
        $der_codon = $3;
        $prob = $4;
        
        $type1 = '';
        $type2 = '';
        
        if (($anc_codon eq 'TTT')||($anc_codon eq 'TTC')) {$anc_AA = 'PHE';}
        if (($anc_codon eq 'TTA')||($anc_codon eq 'TTG')) {$anc_AA = 'LEUa';}
        if (($anc_codon eq 'CTT')||($anc_codon eq 'CTC')||($anc_codon eq 'CTA')||($anc_codon eq 'CTG')) {$anc_AA = 'LEUb';}
        if (($anc_codon eq 'ATT')||($anc_codon eq 'ATC')||($anc_codon eq 'ATA')) {$anc_AA = 'ILE';}
        if (($anc_codon eq 'ATG')) {$anc_AA = 'MET';}
        if (($anc_codon eq 'GTT')||($anc_codon eq 'GTC')||($anc_codon eq 'GTA')||($anc_codon eq 'GTG')) {$anc_AA = 'VAL';}
        if (($anc_codon eq 'TCT')||($anc_codon eq 'TCC')||($anc_codon eq 'TCA')||($anc_codon eq 'TCG')) {$anc_AA = 'SERb';}
        if (($anc_codon eq 'CCT')||($anc_codon eq 'CCC')||($anc_codon eq 'CCA')||($anc_codon eq 'CCG')) {$anc_AA = 'PRO';}
        if (($anc_codon eq 'ACT')||($anc_codon eq 'ACC')||($anc_codon eq 'ACA')||($anc_codon eq 'ACG')) {$anc_AA = 'THR';}
        if (($anc_codon eq 'GCT')||($anc_codon eq 'GCC')||($anc_codon eq 'GCA')||($anc_codon eq 'GCG')) {$anc_AA = 'ALA';}
        if (($anc_codon eq 'TAT')||($anc_codon eq 'TAC')) {$anc_AA = 'TYR';}
        if (($anc_codon eq 'TAA')||($anc_codon eq 'TAG')) {$anc_AA = 'stop';}
        if (($anc_codon eq 'CAT')||($anc_codon eq 'CAC')) {$anc_AA = 'HIS';}
        if (($anc_codon eq 'CAA')||($anc_codon eq 'CAG')) {$anc_AA = 'GLN';}
        if (($anc_codon eq 'AAT')||($anc_codon eq 'AAC')) {$anc_AA = 'ASN';}
        if (($anc_codon eq 'AAA')||($anc_codon eq 'AAG')) {$anc_AA = 'LYS';}
        if (($anc_codon eq 'GAT')||($anc_codon eq 'GAC')) {$anc_AA = 'ASP';}
        if (($anc_codon eq 'GAA')||($anc_codon eq 'GAG')) {$anc_AA = 'GLU';}
        if (($anc_codon eq 'TGT')||($anc_codon eq 'TGC')) {$anc_AA = 'CYS';}
        if (($anc_codon eq 'TGA')) {$anc_AA = 'stop';}
        if (($anc_codon eq 'TGG')) {$anc_AA = 'TRP';}
        if (($anc_codon eq 'CGT')||($anc_codon eq 'CGC')||($anc_codon eq 'CGA')||($anc_codon eq 'CGG')) {$anc_AA = 'ARGb';}
        if (($anc_codon eq 'AGT')||($anc_codon eq 'AGC')) {$anc_AA = 'SERa';}
        if (($anc_codon eq 'AGA')||($anc_codon eq 'AGG')) {$anc_AA = 'ARGa';}
        if (($anc_codon eq 'GGT')||($anc_codon eq 'GGC')||($anc_codon eq 'GGA')||($anc_codon eq 'GGG')) {$anc_AA = 'GLY';}
        
        
        if (($der_codon eq 'TTT')||($der_codon eq 'TTC')) {$der_AA = 'PHE';}
        if (($der_codon eq 'TTA')||($der_codon eq 'TTG')) {$der_AA = 'LEUa';}
        if (($der_codon eq 'CTT')||($der_codon eq 'CTC')||($der_codon eq 'CTA')||($der_codon eq 'CTG')) {$der_AA = 'LEUb';}
        if (($der_codon eq 'ATT')||($der_codon eq 'ATC')||($der_codon eq 'ATA')) {$der_AA = 'ILE';}
        if (($der_codon eq 'ATG')) {$der_AA = 'MET';}
        if (($der_codon eq 'GTT')||($der_codon eq 'GTC')||($der_codon eq 'GTA')||($der_codon eq 'GTG')) {$der_AA = 'VAL';}
        if (($der_codon eq 'TCT')||($der_codon eq 'TCC')||($der_codon eq 'TCA')||($der_codon eq 'TCG')) {$der_AA = 'SERb';}
        if (($der_codon eq 'CCT')||($der_codon eq 'CCC')||($der_codon eq 'CCA')||($der_codon eq 'CCG')) {$der_AA = 'PRO';}
        if (($der_codon eq 'ACT')||($der_codon eq 'ACC')||($der_codon eq 'ACA')||($der_codon eq 'ACG')) {$der_AA = 'THR';}
        if (($der_codon eq 'GCT')||($der_codon eq 'GCC')||($der_codon eq 'GCA')||($der_codon eq 'GCG')) {$der_AA = 'ALA';}
        if (($der_codon eq 'TAT')||($der_codon eq 'TAC')) {$der_AA = 'TYR';}
        if (($der_codon eq 'TAA')||($der_codon eq 'TAG')) {$der_AA = 'stop';}
        if (($der_codon eq 'CAT')||($der_codon eq 'CAC')) {$der_AA = 'HIS';}
        if (($der_codon eq 'CAA')||($der_codon eq 'CAG')) {$der_AA = 'GLN';}
        if (($der_codon eq 'AAT')||($der_codon eq 'AAC')) {$der_AA = 'ASN';}
        if (($der_codon eq 'AAA')||($der_codon eq 'AAG')) {$der_AA = 'LYS';}
        if (($der_codon eq 'GAT')||($der_codon eq 'GAC')) {$der_AA = 'ASP';}
        if (($der_codon eq 'GAA')||($der_codon eq 'GAG')) {$der_AA = 'GLU';}
        if (($der_codon eq 'TGT')||($der_codon eq 'TGC')) {$der_AA = 'CYS';}
        if (($der_codon eq 'TGA')) {$der_AA = 'stop';}
        if (($der_codon eq 'TGG')) {$der_AA = 'TRP';}
        if (($der_codon eq 'CGT')||($der_codon eq 'CGC')||($der_codon eq 'CGA')||($der_codon eq 'CGG')) {$der_AA = 'ARGb';}
        if (($der_codon eq 'AGT')||($der_codon eq 'AGC')) {$der_AA = 'SERa';}
        if (($der_codon eq 'AGA')||($der_codon eq 'AGG')) {$der_AA = 'ARGa';}
        if (($der_codon eq 'GGT')||($der_codon eq 'GGC')||($der_codon eq 'GGA')||($der_codon eq 'GGG')) {$der_AA = 'GLY';}
        
        if ($anc_AA ne $der_AA) {
            if (($anc_AA =~ /LEU/)&&($der_AA =~ /LEU/)) {
                $type1 = 'S';
                $type2 = 'non 3rd pos change';
            }
            elsif (($anc_AA =~ /SER/)&&($der_AA =~ /SER/)) {
                $type1 = 'S';
                $type2 = 'non 3rd pos change';
            }
            elsif (($anc_AA =~ /ARG/)&&($der_AA =~ /ARG/)) {
                $type1 = 'S';
                $type2 = 'non 3rd pos change';
            }
            else {
                $type1 = 'N';
                $type2 = '';
            }
        }
        if ($anc_AA eq $der_AA) {
            $type1 = 'S';
            if ($anc_codon =~ /\w\wT/) {
                if ($der_codon =~ /\w\wC/) {
                    $type2 = 'T->C at 3rd pos';
                }
                if ($der_codon =~ /\w\wA/) {
                    $type2 = 'T->A at 3rd pos';
                }
                if ($der_codon =~ /\w\wG/) {
                    $type2 = 'T->G at 3rd pos';
                }
            }
            if ($anc_codon =~ /\w\wC/) {
                if ($der_codon =~ /\w\wT/) {
                    $type2 = 'C->T at 3rd pos';
                }
                if ($der_codon =~ /\w\wA/) {
                    $type2 = 'C->A at 3rd pos';
                }
                if ($der_codon =~ /\w\wG/) {
                    $type2 = 'C->G at 3rd pos';
                }
            }
            if ($anc_codon =~ /\w\wA/) {
                if ($der_codon =~ /\w\wT/) {
                    $type2 = 'A->T at 3rd pos';
                }
                if ($der_codon =~ /\w\wC/) {
                    $type2 = 'A->C at 3rd pos';
                }
                if ($der_codon =~ /\w\wG/) {
                    $type2 = 'A->G at 3rd pos';
                }
            }
            if ($anc_codon =~ /\w\wG/) {
                if ($der_codon =~ /\w\wT/) {
                    $type2 = 'G->T at 3rd pos';
                }
                if ($der_codon =~ /\w\wC/) {
                    $type2 = 'G->C at 3rd pos';
                }
                if ($der_codon =~ /\w\wA/) {
                    $type2 = 'G->A at 3rd pos';
                }
            }
        }
        
        open (OUT1, ">>list_of_codon_change_with_probability_S_N_Dsim_lineage.txt");
        print (OUT1 "$codon_pos\t$anc_codon\t$der_codon\t$anc_AA\t$der_AA\t$prob\t$type1\t$type2\n");
        close (OUT1);
    }
}

