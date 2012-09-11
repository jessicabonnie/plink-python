=== plink_python ===
Contributors: Jessica Bonnie (jkb4y)
Requires: Python 2.6, PLINK, locuszoom
Version: 1.1.0
Last Updated: 1/17/2012

Set of .py files to perform conditional analysis on data using PLINK.

=== Description ===
plink_python consists of a set of 8 Python program files:

plink_conditional.py
region_wrap.py
parse_log.py
meta_yank.py
assoc_yank.py
map_adapt.py
pc_workhorse.py
pc_toolbox.py

=== BASIC USAGE INFORMATION ===

plink_conditional.py

	=== DESCRIPTION ===
	This program takes a PLINK genotype file, runs it through PLINK,
	illustrates the results via locuszoom, identifies the SNP with the
	most significant (i.e. lowest) p-value, adds this SNP to a condition
	list, and then repeats the process until there are no SNPs with p-values
	below a given threshold of significance.
	Required commands: 	--script <script path>
				--outfolder <folder path>
				--test <test name>
				--chromosome <number>
				--from-mb <start in Mb> --to-mb <end in Mb>
	 		OR	--from-kb <start in kb> --to-kb <end in kb>
			OR	--from-bp <start in bp> --to-bp <end in bp>
				
	Optional commands:	--chrband <chr band>
				--refgene <gene name>
				--flag <word>
				--pbound <number>
				--loop <number>
				--condition-list <file path>
				--pheno <file path>
	
	=== INSTRUCTIONS ===
	Argument and default information for this program can be accessed
	with a call to the program using --help or -h.
	
	== REQUIRED COMMANDS ===
	
	SCRIPT:
	--script <script path>		OR	-s <script path>
	Commands to the program MUST include the path location of a
	script containing the following constructions directly from PLINK:
		--bfile <full path name of binary data files to be fed
		  	into PLINK (without extension)>
		--covar <full path name of covariate file (include extension)>
		Optionally, one can also include the following flags:
		--noweb
		--sex
		--hide-covar
	Scripts CANNOT contain:
		ANY output file information, constructed in PLINK using:
			--out
		ANY type of regional information, constructed in PLINK using:
			--from...
			--to...
			--chromosome
			--window
			--gene
		Association test commands (those will come later):
			For examples, see:
			http://pngu.mgh.harvard.edu/~purcell/plink/anal/ahtml#cc
		Condition List information, constructed in PLINK using:
			--condition-list
		Phenotype commands, constructed in PLINK using
			--pheno
			--all-pheno
			--pheno-name
		
	OUTPUT FOLDER:
	--outfolder <folder path>	OR	-o <folder path>
	Commands to the program MUST include the path location of a Big Mama
	results folder. Within this folder, the program will constuct
	additional folders to separate results by chromosome.
	
	TEST:
	--test <plink test>	OR	-t <plink test>
	Commands to the program MUST include an indication of which PLINK test
	 should be run. Current program configuration will accept the following:
	 	--test logistic
	 	--test linear
	 	--test assoc
	 	--test fisher
	 	--test model
	 	--test mh
	 
	 CHROMOSOME:
	 --chromosome <number>	OR	-c <number>
	 Commands to the program MUST include the number of the chromosome on
	 which the analysis should be run.
	 
	 REGIONAL BOUNDARIES:
	 --from-mb <start in Mb> --to-mb <end in Mb>	OR
	 --from-kb <start in kb> --to-kb <end in kb>	OR
	 --from-bp <start in bp> --to-bp <end in bp>
	 Commands to the program MUST include range information of a region
	 on the chromosome. This can be given in mega-basepairs, kilo-basepairs,
	 or basepairs (or a combination thereof).
	 	__________________________________________
	 
	 === OPTIONAL COMMANDS ===
	 Without these commands, the program will use current default values,
	 which can be learned by using the --help command. To change the
	 defaults, see the advanced user section.
	 
	 CHRBAND:
	 --chrband <chr band>
	 This command uses a chromosomal band name instead of range information
	 to name output files.
	 
	 REFGENE:
	 --refgene <gene name>
	 This command uses a reference gene name instead of range information
	 to name output files.
	 
	 PBOUND:
	 --pbound <highest p-value of interest>
	 This command sets the highest p-value that is considered significant.
	 When there are no more SNPs with p-values below this number,
	 the program terminates.
	 
	 MAXLOOPS:
	 --loop <number> 	OR	-l <number>
	 This command sets the maximum number of times the program should loop
	 through PLINK.
	 
	 FLAG:
	 --flag <word>		OR	-f <word>
	 This command adds a word or phrase (NO SPACES) to the titles of all
	 output files produced during the run.
	 
	 CONDITION-LIST:
	 --condition-list <file path>
	 This command takes the path location of a list of SNPs which should
	 be included in the condition list from the very first loop through
	 PLINK. All SNPs in the file will be copied to ~leastP_SNPs.txt, and
	 any additional SNPs the program identifies will be added there.
	 ANY output files produced through this command will contain '~'
	 in the name.
	 
	 PHENO:
	 --pheno <file path>
	 This command takes the path location of a PLINK phenotype file. The
	 program will perform separate loops for each of the phenotypes in
	 the pheno-file.
	 
	 
	 === EXAMPLES ===
	 Example Script Text:
	 --bfile /home/mst3k/data/soandso
	 --covar /home/mst3k/data/soandso.cov 
	 --hide-covar --sex --noweb
	 
	 Example Commandline Call:
	 python plink_conditional.py --script /home/mst3k/data/soandso.txt
	 --outfolder /home/mst3k/results/SoAndSo/ --test logistic --chromosome 11
	 --refgene INS --from-mb 2.1 --to-mb 2.3
		_______________________________________________

region_wrap.py
	 	 	
	=== DESCRIPTION ===
	This program reads from a file containing regional information and
	performs plink_conditional.py on each of the regions.
	
	Required Commands:	--region-list <file path>
				+ ALL REQUIRED COMMANDS FOR PLINK_CONDITIONAL.PY
				
	Optional Commands:	--cfolder <folder path>
	
	=== INSTRUCTIONS ===
	Argument and default information for this program can be accessed
	with a call to the program using --help or -h.
	
	=== REQUIRED COMMANDS ===
	
	REGION-LIST:
	--region-list <file path>
	Commands to this program MUST include the path to a file containing
	regional information for the regions on which plink_conditional
	should be run. The column headings of the file must include: "gene_chr",
	"region_start", "region_end", and "gene_symbol". The positional 
	information MUST be given in megabase-pairs (Mb).
	
	 
	 PLINK_CONDITIONAL COMMANDS:
	 Commands to this program MUST include any commands that would normally
	 be given to plink_conditional.
	 
	 === OPTIONAL COMMANDS ===
	 
	 CONDITIONAL LIST FOLDER:
	 --cfolder <folder path>
	 This command takes the path location of a folder of condition-lists
	 of SNPs, for naming purposes these should be produced by either
	 meta_yank.py or assoc_yank.py. Each region (e.g. chromosomal band
	 6q72.1 or region JKB on chromosome 6) in the region-list must have
	 a corresponding condition-list in the folder (i.e. 6q72.1.txt or
	 Chr6_JKB.txt). The SNPs in each condition-list will be copied to
	 ~leastP_SNPs.txt and included in the condition-list from the very first
	 loop through PLINK for the associated region. Any additional SNPs the
	 program identifies will be added to the second location
	 (~leastP_SNPs.txt), leaving the original lists in the condition-list
	 folder unchanged. ANY output files produced through this command will
	 contain '~' in the name.
		_____________________________________________________

parse_log.py

	=== DESCRIPTION ===
	This program reads the log files produced during plink_conditional.py
	and creates a summary file containing key information from the logs.
	It also produces a second text file likely of interest only to the
	programmer containing information about the first run through PLINK
	from each log.
	
	Required commands: 	--logfolder <folder path>
						--map <path to ORIGINAL map file>
						--freq <file path>
	
	Optional commands:	--summary <filename NO EXTENSION>
						--runinfo <filename NO EXTENSION>
	
	=== INSTRUCTIONS ===
	Argument and default information for this program can be accessed
	with a call to the program using --help or -h.
	
	=== REQUIRED COMMANDS ===
	
	LOG FOLDER:
	--logfolder <folder path>
	Commands to this program MUST include the path to the logs folder
	which plink_conditional.py created within your Big Mama results folder.
	So, the address of this folder will be Big_Mama's_Path/logs.
	
	MAP FILE:
	--map <file path>
	Commands to this program MUST include the path to the original map file,
	unaltered by map_adapt.
	
	FREQUENCY FILE:
	--freq <file path>
	Commands to this program MUST include the path to the a frequency file for
	the control population produced by plink. For your convenience:
		plink --bfile "data" --filter-controls --freq --make-bed --noweb --out "data"_controls
	
	
	=== OPTIONAL COMMANDS ===
	
	SUMMARY FILE NAME:
	--summary <file name NO EXTENSION>
	The program will write the summary file to the logs folder under the
	name 'log_summary.txt', unless another basename is given.
	
	RUN INFO NAME:
	--runinfo <file name NO EXTENSION>
	The program will write the run info to the logs folder under the name
	'run_info.txt', unless another basename is given.
		_____________________________________________________

meta_yank.py
	=== DESCRIPTION ===
	This program takes a region-list file and a meta-analysis file, both
	with the  specific column heading requirements detailed below.
	It produces a file containing the information of the most significant
	SNP for each of the regions in the region-list. It also produces
	a folder of condition-lists, each containing the most significant
	SNP in a region and named for that region, and a list for each
	region whose lowest p-value is shared by multiple SNPs.
	
	Required Commands:	--outfile <filepath>
						--meta <filepath>
						--region-list <filepath>
	
	=== INSTRUCTIONS ===
	Argument and default information for this program can be accessed
	with a call to the program using --help or -h.
	
	=== REQUIRED COMMANDS ===
	
	OUTPUT FILE:
	--out <file path>
	Commands to this program MUST include the path to the file where the
	results should be written.
	
	META-ANALYSIS FILE:
	--meta <file path>
	Commands to this program MUST include the path to a meta-analysis file
	containing the relevant SNP information. The column headings of the file
	must include: "CHR", "MarkerName", "P-value","POS" (bp-position),
	"Zscore".
	
	REGION-LIST:
	--region-list <file path> or -r <file path>
	Commands to this program MUST include the path to a region-list file. The
	column headings of the file must include: "gene_chr", "region_start",
	"region_end", and "gene_symbol". The positional information MUST be given
	in megabase-pairs (Mb).

		_____________________________________________________
assoc_yank.py
	=== DESCRIPTION ===
	This program takes a region-list file and a PLINK association file, both
	with the  specific column heading requirements detailed below.
	It produces a file containing the information of the most significant
	SNP for each of the regions in the region-list. It also produces
	a folder of condition-lists, each containing the most significant
	SNP in a region (if it could be determined)
	and named for that region, and a list for each region whose lowest
	p-value is shared by multiple SNPs.
	
	Required Commands:	--outfile <filepath>
				--assoc <filepath>
				--region-list <filepath>
				--bp-form <'mb' OR 'kb' OR 'bp'>
	
	=== INSTRUCTIONS ===
	Argument and default information for this program can be accessed
	with a call to the program using --help or -h.
	
	=== REQUIRED COMMANDS ===
	
	OUTPUT FILE:
	--out <file path>
	Commands to this program MUST include the path to the file where the
	results should be written.
	
	ASSOCIATION FILE:
	--assoc <file path>
	Commands to this program MUST include the path to a PLINK association file
	containing the relevant SNP information. The column headings of the file
	must include: "CHR", "SNP", "P", "BP" (bp-position).
	
	REGION-LIST:
	Commands to this program MUST include the path to a region-list file. The
	column headings of the file must include: "gene_chr", "region_start",
	"region_end", and "gene_symbol".
	
	POSITIONAL UNITS:
	--bp-form <'mb' OR 'kb' OR 'bp'>
	Commands to this program MUST include one of three flags indicating
	what unit is used to list positional information in the region list.

	
		_____________________________________________________
	
map_adapt.py
	=== DESCRIPTION ===
	This program takes a map file (e.g.'.bim') and produces a copy of
	the file replacing any names of non-duplicate SNPs that are not
	in 'rs' form with chr<chromosome>:<basepair position> names. The 
	original map file is saved to another name, while the new one
	assumes the name of the original. The program writes a list of the
	duplicate SNPs to a new file in the same folder at the map file.
	
	Required commands:	--map <file path>
	
	=== INSTRUCTIONS ===
	Argument information for this program can be accessed with a call
	to the program using --help or -h.
	
	=== REQUIRED COMMANDS ===
	
	MAP FILE:
	--map <path to map file>
	Commands to the program MUST include the path to a map file.
		__________________________________________________
	
pc_workhorse.py
	=== DESCRIPTION ===
	This program does all the work that is instantiated by a call to
	plink_conditional.py.
	
	=== INSTRUCTIONS ===
	
	DO NOT MAKE CALLS TO PC_WORKHORSE.PY
	

pc_toolbox.py
	=== DESCRIPTION ===
	This program contains functions which are used by multiple programs
	in the suite.
	
	=== INSTRUCTIONS ===
	
	DO NOT MAKE CALLS TO PC_TOOLBOX.PY


=== ADVANCED USAGE INFORMATION ===

	ONLY Users who know what they are doing:
		CHANGING DEFAULTS
		It might be easier for you to change some of the defaults at
		the top	of the program file than type in all of the command
		flags every time. If you would like, you may do this for
		plink_conditional.py or parse_log.py.
		DO NOT CHANGE DEFAULTS IN PC_WORKHORSE.PY.
		Open the program file in a text editor. The default values are
		set beneath the list of globals. The ones in which you may be
		interested are marked between rows of asterisks (*). DO NOT
		alter anything (i.e. the name of the default) on the left side
		of the equals sign, only change the value on the right side.
	PLINK Constructions:
		plink_conditional.py only cares about 2 things:
		  1) The column names of the PLINK results file. The program
		     will look for columns titled 'P' and 'SNP' in order to
		     perform its search for the most significant SNP and then
		     use those column titles to instruct locus zoom. If a PLINK
		     construction results in differently titled columns, there
		     will be problems.
		  2) The file extension of the PLINK results file containing the
		     SNP and P-value information. If a PLINK construction alters
		     the file extension, there will be problems.
	



=== CHANGE LOG ===

**10/26/2011:
Changes made to condition-list functionality to avoid overwriting output produced during a
"traditional" run of plink_conditional without providing a condition-list.

Bug Fixed: condition-list SNPs accidentally overwritten in leastP_SNPs file after first run.

**10/27/2011:
Changes made to logging functionality to permit condition_list existance and flag information
to be extracted during parse_log.

**10/28/2011
Changed name of worker program from plink_association.py to pc_workhorse.py to better deter
user from opening it.

**10/30/2011
Updated Advanced Usage section of README file to further explain PLINK construction options.

**11/1/2011 - 11/13/2011
Additional program, region_wrap, added to package in order to facilitate the
reading of regional information by plink_conditional from a region-list.

Changes made to plink_conditional to allow for use of phenotype related PLINK
functionality during tests.

**11/22/2011 - NOW v1.0.3
Bug Fixed: region_wrap misinterpreting command flags.
Bug Fixed: filenames of locuszoom pdfs contained '.assoc'
Additional programs, meta_yank and assoc_yank, added to package to facilitate production of
condition-lists.

**1/17/2012 - NOW v1.1.0
Additional program, pc_toolbox, added to package to improve modularity
Changes made to pc_workhorse to allow for identification of most significant SNP
after final run.
Changes made to entire suite to allow use of chromosomal band information during naming.
Changes made to pc_workhorse to include y-axis label in locuszoom plots
Changes made to meta_yank to include frequency information for selected SNPs
Changes made to parse_log to include frequency information for significant SNPs

**2/2/2012
Changes made to pc_workhorse to allow choice between to SNPs with p=0 using |t-statistic|.
