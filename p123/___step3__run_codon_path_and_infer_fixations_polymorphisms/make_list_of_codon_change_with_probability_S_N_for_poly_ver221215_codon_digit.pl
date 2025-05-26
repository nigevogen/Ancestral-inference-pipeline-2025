open (IN1,"list_of_codon_change_with_probability_Dsim_poly_codon_digit.txt");
@str1 = <IN1>;

open (OUT1, ">>list_of_codon_change_with_probability_S_N_Dsim_poly_codon_digit.txt");

foreach $_ (@str1) {
    if ($_ =~ /(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d.*)\n/) {
        $codon_pos = $1;
        $anc_codon = $2;
        $der_codon = $3;
        $freq = $4;
        $prob = $5;
        
        $type1 = '';
        $type2 = '';
        
        if (($anc_codon == 1)||($anc_codon == 2)) {$anc_AA = 'PHE';}
        if (($anc_codon == 3)||($anc_codon == 4)) {$anc_AA = 'LEU2f';}
        if (($anc_codon == 5)||($anc_codon == 6)||($anc_codon == 7)||($anc_codon == 8)) {$anc_AA = 'LEU4f';}
        if (($anc_codon == 9)||($anc_codon == 10)||($anc_codon == 11)) {$anc_AA = 'ILE';}
        if (($anc_codon == 12)) {$anc_AA = 'MET';}
        if (($anc_codon == 13)||($anc_codon == 14)||($anc_codon == 15)||($anc_codon == 16)) {$anc_AA = 'VAL';}
        if (($anc_codon == 17)||($anc_codon == 18)||($anc_codon == 19)||($anc_codon == 20)) {$anc_AA = 'SER4f';}
        if (($anc_codon == 21)||($anc_codon == 22)||($anc_codon == 23)||($anc_codon == 24)) {$anc_AA = 'PRO';}
        if (($anc_codon == 25)||($anc_codon == 26)||($anc_codon == 27)||($anc_codon == 28)) {$anc_AA = 'THR';}
        if (($anc_codon == 29)||($anc_codon == 30)||($anc_codon == 31)||($anc_codon == 32)) {$anc_AA = 'ALA';}
        if (($anc_codon == 33)||($anc_codon == 34)) {$anc_AA = 'TYR';}
        if (($anc_codon == 35)||($anc_codon == 36)) {$anc_AA = 'stop';}
        if (($anc_codon == 37)||($anc_codon == 38)) {$anc_AA = 'HIS';}
        if (($anc_codon == 39)||($anc_codon == 40)) {$anc_AA = 'GLN';}
        if (($anc_codon == 41)||($anc_codon == 42)) {$anc_AA = 'ASN';}
        if (($anc_codon == 43)||($anc_codon == 44)) {$anc_AA = 'LYS';}
        if (($anc_codon == 45)||($anc_codon == 46)) {$anc_AA = 'ASP';}
        if (($anc_codon == 47)||($anc_codon == 48)) {$anc_AA = 'GLU';}
        if (($anc_codon == 49)||($anc_codon == 50)) {$anc_AA = 'CYS';}
        if (($anc_codon == 51)) {$anc_AA = 'stop';}
        if (($anc_codon == 52)) {$anc_AA = 'TRP';}
        if (($anc_codon == 53)||($anc_codon == 54)||($anc_codon == 55)||($anc_codon == 56)) {$anc_AA = 'ARG4f';}
        if (($anc_codon == 57)||($anc_codon == 58)) {$anc_AA = 'SER2f';}
        if (($anc_codon == 59)||($anc_codon == 60)) {$anc_AA = 'ARG2f';}
        if (($anc_codon == 61)||($anc_codon == 62)||($anc_codon == 63)||($anc_codon == 64)) {$anc_AA = 'GLY';}
        
        
        if (($der_codon == 1)||($der_codon == 2)) {$der_AA = 'PHE';}
        if (($der_codon == 3)||($der_codon == 4)) {$der_AA = 'LEU2f';}
        if (($der_codon == 5)||($der_codon == 6)||($der_codon == 7)||($der_codon == 8)) {$der_AA = 'LEU4f';}
        if (($der_codon == 9)||($der_codon == 10)||($der_codon == 11)) {$der_AA = 'ILE';}
        if (($der_codon == 12)) {$der_AA = 'MET';}
        if (($der_codon == 13)||($der_codon == 14)||($der_codon == 15)||($der_codon == 16)) {$der_AA = 'VAL';}
        if (($der_codon == 17)||($der_codon == 18)||($der_codon == 19)||($der_codon == 20)) {$der_AA = 'SER4f';}
        if (($der_codon == 21)||($der_codon == 22)||($der_codon == 23)||($der_codon == 24)) {$der_AA = 'PRO';}
        if (($der_codon == 25)||($der_codon == 26)||($der_codon == 27)||($der_codon == 28)) {$der_AA = 'THR';}
        if (($der_codon == 29)||($der_codon == 30)||($der_codon == 31)||($der_codon == 32)) {$der_AA = 'ALA';}
        if (($der_codon == 33)||($der_codon == 34)) {$der_AA = 'TYR';}
        if (($der_codon == 35)||($der_codon == 36)) {$der_AA = 'stop';}
        if (($der_codon == 37)||($der_codon == 38)) {$der_AA = 'HIS';}
        if (($der_codon == 39)||($der_codon == 40)) {$der_AA = 'GLN';}
        if (($der_codon == 41)||($der_codon == 42)) {$der_AA = 'ASN';}
        if (($der_codon == 43)||($der_codon == 44)) {$der_AA = 'LYS';}
        if (($der_codon == 45)||($der_codon == 46)) {$der_AA = 'ASP';}
        if (($der_codon == 47)||($der_codon == 48)) {$der_AA = 'GLU';}
        if (($der_codon == 49)||($der_codon == 50)) {$der_AA = 'CYS';}
        if (($der_codon == 51)) {$der_AA = 'stop';}
        if (($der_codon == 52)) {$der_AA = 'TRP';}
        if (($der_codon == 53)||($der_codon == 54)||($der_codon == 55)||($der_codon == 56)) {$der_AA = 'ARG4f';}
        if (($der_codon == 57)||($der_codon == 58)) {$der_AA = 'SER2f';}
        if (($der_codon == 59)||($der_codon == 60)) {$der_AA = 'ARG2f';}
        if (($der_codon == 61)||($der_codon == 62)||($der_codon == 63)||($der_codon == 64)) {$der_AA = 'GLY';}
        
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
        print (OUT1 "$codon_pos\t$anc_codon\t$der_codon\t$anc_AA\t$der_AA\t$freq\t$prob\t$type1\t$type2\n");
    }
}
close (OUT1);

open (IN1,"list_of_codon_change_with_probability_Dmel_poly_codon_digit.txt");
@str1 = <IN1>;

open (OUT1, ">>list_of_codon_change_with_probability_S_N_Dmel_poly_codon_digit.txt");

foreach $_ (@str1) {
    if ($_ =~ /(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d.*)\n/) {
        $codon_pos = $1;
        $anc_codon = $2;
        $der_codon = $3;
        $freq = $4;
        $prob = $5;
        
        $type1 = '';
        $type2 = '';
        
        if (($anc_codon == 1)||($anc_codon == 2)) {$anc_AA = 'PHE';}
        if (($anc_codon == 3)||($anc_codon == 4)) {$anc_AA = 'LEU2f';}
        if (($anc_codon == 5)||($anc_codon == 6)||($anc_codon == 7)||($anc_codon == 8)) {$anc_AA = 'LEU4f';}
        if (($anc_codon == 9)||($anc_codon == 10)||($anc_codon == 11)) {$anc_AA = 'ILE';}
        if (($anc_codon == 12)) {$anc_AA = 'MET';}
        if (($anc_codon == 13)||($anc_codon == 14)||($anc_codon == 15)||($anc_codon == 16)) {$anc_AA = 'VAL';}
        if (($anc_codon == 17)||($anc_codon == 18)||($anc_codon == 19)||($anc_codon == 20)) {$anc_AA = 'SER4f';}
        if (($anc_codon == 21)||($anc_codon == 22)||($anc_codon == 23)||($anc_codon == 24)) {$anc_AA = 'PRO';}
        if (($anc_codon == 25)||($anc_codon == 26)||($anc_codon == 27)||($anc_codon == 28)) {$anc_AA = 'THR';}
        if (($anc_codon == 29)||($anc_codon == 30)||($anc_codon == 31)||($anc_codon == 32)) {$anc_AA = 'ALA';}
        if (($anc_codon == 33)||($anc_codon == 34)) {$anc_AA = 'TYR';}
        if (($anc_codon == 35)||($anc_codon == 36)) {$anc_AA = 'stop';}
        if (($anc_codon == 37)||($anc_codon == 38)) {$anc_AA = 'HIS';}
        if (($anc_codon == 39)||($anc_codon == 40)) {$anc_AA = 'GLN';}
        if (($anc_codon == 41)||($anc_codon == 42)) {$anc_AA = 'ASN';}
        if (($anc_codon == 43)||($anc_codon == 44)) {$anc_AA = 'LYS';}
        if (($anc_codon == 45)||($anc_codon == 46)) {$anc_AA = 'ASP';}
        if (($anc_codon == 47)||($anc_codon == 48)) {$anc_AA = 'GLU';}
        if (($anc_codon == 49)||($anc_codon == 50)) {$anc_AA = 'CYS';}
        if (($anc_codon == 51)) {$anc_AA = 'stop';}
        if (($anc_codon == 52)) {$anc_AA = 'TRP';}
        if (($anc_codon == 53)||($anc_codon == 54)||($anc_codon == 55)||($anc_codon == 56)) {$anc_AA = 'ARG4f';}
        if (($anc_codon == 57)||($anc_codon == 58)) {$anc_AA = 'SER2f';}
        if (($anc_codon == 59)||($anc_codon == 60)) {$anc_AA = 'ARG2f';}
        if (($anc_codon == 61)||($anc_codon == 62)||($anc_codon == 63)||($anc_codon == 64)) {$anc_AA = 'GLY';}
        
        
        if (($der_codon == 1)||($der_codon == 2)) {$der_AA = 'PHE';}
        if (($der_codon == 3)||($der_codon == 4)) {$der_AA = 'LEU2f';}
        if (($der_codon == 5)||($der_codon == 6)||($der_codon == 7)||($der_codon == 8)) {$der_AA = 'LEU4f';}
        if (($der_codon == 9)||($der_codon == 10)||($der_codon == 11)) {$der_AA = 'ILE';}
        if (($der_codon == 12)) {$der_AA = 'MET';}
        if (($der_codon == 13)||($der_codon == 14)||($der_codon == 15)||($der_codon == 16)) {$der_AA = 'VAL';}
        if (($der_codon == 17)||($der_codon == 18)||($der_codon == 19)||($der_codon == 20)) {$der_AA = 'SER4f';}
        if (($der_codon == 21)||($der_codon == 22)||($der_codon == 23)||($der_codon == 24)) {$der_AA = 'PRO';}
        if (($der_codon == 25)||($der_codon == 26)||($der_codon == 27)||($der_codon == 28)) {$der_AA = 'THR';}
        if (($der_codon == 29)||($der_codon == 30)||($der_codon == 31)||($der_codon == 32)) {$der_AA = 'ALA';}
        if (($der_codon == 33)||($der_codon == 34)) {$der_AA = 'TYR';}
        if (($der_codon == 35)||($der_codon == 36)) {$der_AA = 'stop';}
        if (($der_codon == 37)||($der_codon == 38)) {$der_AA = 'HIS';}
        if (($der_codon == 39)||($der_codon == 40)) {$der_AA = 'GLN';}
        if (($der_codon == 41)||($der_codon == 42)) {$der_AA = 'ASN';}
        if (($der_codon == 43)||($der_codon == 44)) {$der_AA = 'LYS';}
        if (($der_codon == 45)||($der_codon == 46)) {$der_AA = 'ASP';}
        if (($der_codon == 47)||($der_codon == 48)) {$der_AA = 'GLU';}
        if (($der_codon == 49)||($der_codon == 50)) {$der_AA = 'CYS';}
        if (($der_codon == 51)) {$der_AA = 'stop';}
        if (($der_codon == 52)) {$der_AA = 'TRP';}
        if (($der_codon == 53)||($der_codon == 54)||($der_codon == 55)||($der_codon == 56)) {$der_AA = 'ARG4f';}
        if (($der_codon == 57)||($der_codon == 58)) {$der_AA = 'SER2f';}
        if (($der_codon == 59)||($der_codon == 60)) {$der_AA = 'ARG2f';}
        if (($der_codon == 61)||($der_codon == 62)||($der_codon == 63)||($der_codon == 64)) {$der_AA = 'GLY';}
        
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
        
        print (OUT1 "$codon_pos\t$anc_codon\t$der_codon\t$anc_AA\t$der_AA\t$freq\t$prob\t$type1\t$type2\n");
    }
}
close (OUT1);
