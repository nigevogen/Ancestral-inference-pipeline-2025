/*--------------------------------------------------*/
process_choice	= 331	/* */
MainFolder	= "../04_output/tmp_result"	/* */
CreateFldOpt	= 1	/* */
bin_gn_dat_list_opt     = 1     /* */
ConcatList	= "../../03_seq_dat/sample_data/sample_concat_list_codons_nucleotides/bin_dat_3_index_lists"	/* */
MfaFld_in		= "../../03_seq_dat/sample_data/sample_sequences"	/* */
MfaFld_out	= "2_alncds_cc.MFA"	/* */
NumErrorCodesOK	= 6	/* */
ErrorCodesOK	= 0,1,2,4,8,16	/* */
rootInputFolder	: "../10.input/"	/* */
rootOutputFolder	: "../20.output"	/* */
dirSep 			: "/"				/* */
/*--------------------------------------------------*/
fileOpenMode		: TEXT				/* */
newLine			: LF				/* */
overWriteOption	: dontOverWrite		/* */
/*--------------------------------------------------*/
