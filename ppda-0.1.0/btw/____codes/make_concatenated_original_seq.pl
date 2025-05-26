#read list of the sequences to be concatenated
open (IN1,"../_seq_list_to_be_concatenated/seq_list");
@str1 = <IN1>;

#search the sequence listed in the list
for ($line=0;$line<=$#str1;$line++) {
    if (@str1[$line] =~/^(.*)/) {
        $ID = $1;
        
        open (IN2,"../_seqs_folder/$ID");
        @str2 = <IN2>;
        
        #separate sequence line and name line
        $count = 1;
        foreach $_ (@str2) {
            if ($_ =~/\>.*\n/) {
                $index1[$count] = $_;
            }
            if (($_ !~/\>/) && ($_ =~/^A||^T||^G||^C||^N||^\-/)) {
                $seq1 = $_;
                chomp($seq1);
                
                #concatenation
                $seq1[$count] = $seq1[$count].$seq1;
                $count++;
            }
        }
    }
}
open (OUT1, ">>concatenated_original_alignment.txt");

for ($n=1;$n<=$count-1;$n++) {
    print (OUT1 "$index1[$n]");
    print (OUT1 "$seq1[$n]\n");
}
print (OUT1 "\n");
close (OUT1);
