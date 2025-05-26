open (IN1,"joint_probability_of_ancestral_codon_configuration_from_BASEML_and_BTW_wo_configuration_with_stop_codons_codon_digit.txt");
@str1 = <IN1>;
open (IN3,"bin_dat_1_index_lists");
@str3= <IN3>;


$st = 0;
$bookmark_str1 = 0;
foreach $_ (@str3) {
    if ($_ =~ /^(Dmel.*)\.mpspye\.aln\n/) {
    $geneID = $1;
        if (-e "m14s21ye_r624_CDS_aln_201217HY_without_N_G_M_stopcodon_filter_by_three_markers_1st_pos_filter_collapse/$geneID\.mpspye\.aln") {
            open (IN2,"m14s21ye_r624_CDS_aln_201217HY_without_N_G_M_stopcodon_filter_by_three_markers_1st_pos_filter_collapse/$geneID\.mpspye\.aln");
            @str2 = <IN2>;
            $check = 0;
            if (@str2[1] =~ /^\w+/) {
                $check = 1;
            }
            
            if ($check == 1) {
                
                open (OUT1, ">>codon_count_at_m_each_gene_codon_digit.txt");
                print (OUT1 "$geneID\t");
                close (OUT1);
                
                open (OUT1, ">>codon_count_at_s_each_gene_codon_digit.txt");
                print (OUT1 "$geneID\t");
                close (OUT1);
                
                for ($n1=1;$n1<=64;$n1++) {
                    ${'count_m_'.$n1} = 0;
                    ${'count_s_'.$n1} = 0;
                }
                
                $length = length @str2[1];
                $length = $length - 1;
                $ed = $st + $length - 1;
                
                for ($line=$bookmark_str1;$line<=$#str1;$line++) {
                    if (@str1[$line] =~ /^(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d.*)\n/) {
                        $codon_pos = $1;
                        $ms = $2;
                        $ty = $3;
                        $m = $4;
                        $s = $5;
                        $prob = $6;
                        
                        if ($codon_pos > $ed) {
                            $st = $st + $length;
                            $bookmark_str1 = $line;
                            last;
                        }
                        if (($codon_pos >= $st)&&($codon_pos <= $ed)) {
                            
                            ${'count_m_'.$ms} = ${'count_m_'.$ms} + $prob;
                            ${'count_s_'.$ms} = ${'count_s_'.$ms} + $prob;
                        }
                    }
                }
                for ($n1=1;$n1<=64;$n1++) {
                    open (OUT1, ">>codon_count_at_m_each_gene_codon_digit.txt");
                    print (OUT1 "${'count_m_'.$n1}\t");
                    close (OUT1);
                }
                for ($n1=1;$n1<=64;$n1++) {
                    open (OUT1, ">>codon_count_at_s_each_gene_codon_digit.txt");
                    print (OUT1 "${'count_s_'.$n1}\t");
                    close (OUT1);
                }
                open (OUT1, ">>codon_count_at_m_each_gene_codon_digit.txt");
                print (OUT1 "\n");
                close (OUT1);
                open (OUT1, ">>codon_count_at_s_each_gene_codon_digit.txt");
                print (OUT1 "\n");
                close (OUT1);
            }
        }
    }
}

