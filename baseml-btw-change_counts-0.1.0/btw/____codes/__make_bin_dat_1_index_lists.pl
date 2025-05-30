#make bin_dat_1_index_list by reading the seq_list file

open (IN1,"../_seq_list_to_be_concatenated/seq_list");
@str1 = <IN1>;

$size = 0;
foreach $_ (@str1) {
    if ($_ =~ /.*\n/) {
        $size++;
    }
}

#this is the required description to run our BASEML pipeline

open (OUT1, ">>../bin_dat_1_index_lists");
printf (OUT1 "1\n");
printf (OUT1 "prefix: \n");
printf (OUT1 "suffix: \n");
printf (OUT1 ">concate_seq\t$size\n");
close (OUT1);

foreach $_ (@str1) {
    
    open (OUT1, ">>../bin_dat_1_index_lists");
    printf (OUT1 "$_");
    close (OUT1);
}

