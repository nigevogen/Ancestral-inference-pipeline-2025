/*-----------------------------------------------------------------------------------------------------------------*/
/* folder and file names and parameters */
/*---------------------------------------------------------------------------------------*/
process_choice				= 335						/* Process shifting option */
MainFolder				= "../04_output/tmp_result"		/* main path of the data folder */
/* options and constant num */
CreateFldOpt			= 0			/* 0:not create a new folder; 1:create a new folder */
CodNumCut				= 0			/* cut short genes */
SmallValCut			    = 0.0000	/* threshold for cutting genes */
SpeciesNum				= 4			/* species number */
SBR_opt				    = 0			/* single best reconstruction: 2, SBR does not use even prob; 1, SBR use even prob; 0, WAP*/
RangeStart				= -9		/* read from the RangeStart-th gene  */
RangeEnd				= -9		/* read to the RangeEnd-th gene */
FormatOpt				= 2			/* option used to select different format: 0, HM; 1, HM01; 2, HM01s. */
debug					= 0			/* debug option: 0, not debug; 1, detail; 2, debug with halts. */
bin_gn_dat_list_opt     = 1   		/* bin gene dat list output option (with codon num per gene) */
/* Codon Adjust Option */
StopCodon_opt			= 3		/* 0:no change; 1:cut terminal cod; 2:pro=0; 3:redistribute prob */
mHits_opt				= 3		/* 0:no change; 1:cut terminal cod; 2:pro=0; 3:redistribute prob */
/*---------------------------------------------------------------------------------------*/

MainRstBMLFld		= "5_alncds_filter_rst"		/* CDS (DNA) 3 position baseml output folder */
RstBMLPos1			= "pos1_fltrst"						/* CDS (DNA) 1st position baseml reconstruction file */
RstBMLPos2			= "pos2_fltrst"						/* CDS (DNA) 2nd position baseml reconstruction file */
RstBMLPos3			= "pos3_fltrst"						/* CDS (DNA) 3rd position baseml reconstruction file */
BinGnDatList		= "0.bin_gene_dat_list"	/* bin gene dat list file */

MainBML45CodBrnchSummFld	= "8_alncds_HM"	/* CDS (DNA) 3 position codon branch prob. summary */
BMLDSBrnchHM				= "alncds_model1_HM"	/* CDS (DNA) 3 position codon change prob. summary, HitMetrix */

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