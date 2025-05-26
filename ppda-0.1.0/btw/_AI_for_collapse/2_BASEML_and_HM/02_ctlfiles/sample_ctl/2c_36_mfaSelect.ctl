/*-----------------------------------------------------------------------------------------------------------------*/
/* folder and file names and parameters */
/*---------------------------------------------------------------------------------------*/
process_choice				= 332			/* 1=codon pref test, 2=tissue MCUsp aa props, 3=tissue MCU v expr rs, 5=plot MCUsp in each tissue, see main() for others */
MainFolder				= "../04_output/tmp_result"		/* mod_mss; main path of the data folder */
/* options and constant num */
CreateFldOpt			= 1			/* 0:not create a new folder; 1:create a new folder */
RangeStart				= 1		/* output from the RangeStart-th gene  */
RangeEnd				= 30		/* output to the RangeEnd-th gene */
/*---------------------------------------------------------------------------------------*/

MfaFld				= "2_alncds_cc.MFA"		/* mfa folder (input) */
FltMfaFld			= "2c_alncds_cc_slct.MFA"  /* mfa folder (output) */

/*---------------------------------------------------------------------------------------*/
sample_long					= 1													/* brief explanation */
sample_string				= "../../../00_data/dmel/528r_X/gene_map_HA2.txt"	/* brief explanation */

NumErrorCodesOK				= 6					/* The number of errorcodes that will be allowed in the final ffn file */
ErrorCodesOK				= 0,1,2,4,8,16   	/* The error codes that will be pulled through to the output (024 is standard) */
												/* If it is desired that only genes with no errors are outputted */
												/*	NumErrorCodesOK =1 and ErrorCodesOk = 0			*/
												/*  0 : No error was present in the gene			*/
												/*  1 : Basenum was a non-multiple of 3 ***			*/
												/*  2 : Non ATG start codon							*/
												/*  4 : Bad termination codon						*/
												/*  8 : Non ATCG presence ***						*/
												/*  16: Premature Stop codon presence ***			*/	

/*---------------------------------------------------------------------------------------*/
/*these values usually don't have to be changed*/
/*---------------------------------------------------------------------------------------*/
rootInputFolder		: "../10.input/"					/* path+name of the root input folder */
rootOutputFolder	: "../20.output/"				/* path+name of the root output folder /Users/boyang/Documents/01_HA_Xcode/99_BL_ctl/20.output/ */
dirSep 				: "/"							/* the directory separator to be used */
/*---------------------------------------------------------------------------------------*/
/*new file preferences         Not being currently used*/
/*---------------------------------------------------------------------------------------*/
fileOpenMode		: TEXT						/*the new file open mode BINARY or TEXT*/
newLine				: LF						/*the new line encode CR, LF or CRLF where CR = \r LF = \n*/
overWriteOption		: dontOverWrite				/*overWrite = files will be overwritten if already exist
												  dontOverWrite = program will exit if files already exist
												  askUser = program will ask user as to what should be done*/
/*---------------------------------------------------------------------------------------*/
/*donot change the names of the parameters as the program will be checking them to
match parameters. order of the parameters doesnt matter. both styles of c99 commenting
are supported and comments can run into multiple lines. The fields with a colon ":" as the
separator indicates values that are less frequently changed and the fields with a "=" are
values that would probably need to be changed for every new run. FileNames, foldernames and
paths are enclosed in quotes */