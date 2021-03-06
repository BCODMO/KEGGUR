# KEGGUR
A suite of perl scripts for generating "good enough" phylogenetic databases using the Kyoto Encyclopedia of Genes and Genomes (KEGG)

An archival copy of this code repository, and metadata can be found at BCO-DMO: https://www.bco-dmo.org/dataset/636508

NOTE: view this file as raw text, not html, as syntax assignments are represented using <> notation.

KEGGUR v 0.1

readme.txt

12 August 2017

OVERVIEW:

KEGGUR is a perl-based software package designed to assist in detecting and analyzing protein-coding genes in metagenomic and metatranscriptomic datasets. It consists of 3 stand-alone main scripts as well as a number of specialized scripts specifically designed for the pipeline used in the Morris et al. 2015 J Plankton Research paper.

MAIN PROGRAMS:
getorgs.pl
keggur.pl
manalign.pl

SPECIALIZED PROGRAMS:
fixnames.pl
makegroups.pl
metaprocess.pl
blastseqs.postblast.pl


getorgs.pl

Description: This script retrieves a current list of all reference genomes contained in the Kyoto Encyclopedia of Genes and Genomes (KEGG) and outputs four files.  First is a simple list of all the 3- or 4-character organism codes (kegg.orgs).  Second is a key file with those codes followed by the full name for each organism (kegg.names).  Third is a master taxonomy file containing the organism code followed by the entire taxonomic assignment from Domain to Species for each organism (kegg.tax). Fourth is a logfile containing the keggur version number, date of operation, and various outputs (getorgs.logfile)

Dependencies: getorgs.pl requires WGET to be available in the same path that getorgs.pl is called from.  Also requires POSIX perl module.

Syntax: getorgs.pl is run without arguments.


keggur.pl

Description: This script provides keggur's main functionality.  Given an orthology code for a functional gene, keggur.pl gets a list of all instances of that gene in the keggur database.  Each instance is associated with an organism (found in the output taxonomy file from getorgs.pl).  Keggur downloads the nucleotide sequences of all of these genes and eliminates any identical sequences, then translates them all to amino acids and aligns them using MUSCLE.  The MUSCLE alignment program is run multiple times until it either converges on a state that doesn't change between iterations or converges on a loop of equally good alignments.  In the latter case, keggur selects the shortest alignment as the preferable one.  After aligning the sequences as amino acids, keggur converts the alignment back to nucleotides by reference to the original nucleotide sequences.  Next, keggur produces a conservative taxonomy file (in MOTHUR format) that makes sure the taxonomies of any organisms containing identical sequences are also identical and, if they don't, eliminates taxonomic assignments below the level of disagreement.  Finally, keggur generates a BLAST database for the functional gene.

Keggur can be given a maximum number of muscle iterations to use (if the alignments are taking too long, for instance).  If keggur doesn't decide on an ideal alignment in this number of iterations, it outputs the last alignment (which is probably still pretty good).

Also, the processing time can be reduced by telling keggur to use OTU clustering on the sequences.  In this case, keggur is provided a cutoff similarity and lumps sufficiently similar sequences into a single OTU, then selects a representative sequence to output in the final files.  Note that keggur's taxonomy file will reflect the results of OTU clustering and will give a good-faith representation of the taxa represented by each representative sequence.
	
Dependencies: keggur.pl requires WGET, MOTHUR, MUSCLE, and BLAST+ in the same path from which keggur.pl is called.  It uses the POSIX, Getopt::Std, List::Util, and File::Path perl modules.

Syntax:
	perl keggur.pl <orthology code> <file/gene name> <output directory> <maximum number of muscle iterations to use> <master taxonomy file from getorgs.pl> <whether or not to use OTU clustering> <cutoff identity for clustering, between 0 and 1>
	
	
manalign.pl

Description: This script (MANual ALIGNment) is for inserting short-read BLAST hits identified in metagenomes or metatranscriptomes using the keggur-generated databases.  Rather than use an alignment algorithm, each blast hit is placed into the alignment by direct comparison to its best reference match identified by BLAST.  Where the reference has a gap, the aligned sequence will also have a gap.  In the case that a new gap is placed in the reference sequence by BLAST, manalign.pl will insert a gap at that position in every sequence in the reference alignment.

Note that manalign.pl expects to be given a blastx output file in a very specific format.  When blastx is called, it should include the tag:

-outfmt "10 qseqid sseqid qseq sseq evalue bitscore length gaps pident qcovs qstart qend sstart send qframe mismatch"

Also manalign.pl expects to be given amino acid data, but outputs nucleotide alignments.

Note that manalign.pl is designed to work with the output file generated by blastseqs.postblast.pl (see below).

Dependencies: manalign.pl requires an output .csv file generated by blastx.  It also uses the POSIX perl module.

Syntax:
	perl manalign.pl <KEGGUR reference nucleotide alignment> <fasta file containing nucleotide sequences to be inserted into the reference alignment> <BLAST csv file containing each sequence in the previously mentioned fasta file> <name of output file>


fixnames.pl

Description:  This script finds KEGG organism codes in a file by looking between two delimiting strings.  It then replaces those codes with the organism's name and, if a KEGGUR taxonomy file is given, will also insert the organism's entire taxonomy.  We used this file to replace codes in NEWICK phylogenetic tree files with taxonomic data to facilitate "eyeball" analysis of the trees.

Dependencies: none

Syntax:
	perl fixnames.pl <file containing codes to be changed> <output file> <file containing keg organism codes -- probably kegg.names from getorgs.pl> <delimiter BEFORE the name code> <delimiter AFTER the name code> <optional taxonomy file -- probably kegg.tax from getorgs.pl>
	

makegroups.pl

Description: This is a simple script to extract groups information from a master fasta file containing multiple samples and output a mothur groups file, which can then be used either to break the master file into individual sample file using mothur.  This was necessary for us because the CAMERA metatranscriptome files compiled all samples from a single project into a single large fasta file.  The group name is discovered using a regular expression provided when makegroups.pl is called that extracts the group a sequence is associated with from the sequence's fasta name line.

Dependencies:  This script requires the POSIX, Getopt::Std, and File::Path perl modules.

Syntax:
	perl makegroups.pl <sequence filename to be analyzed> <regular expression to be used> <output filename> <logfile name>
	
	
metaprocess.pl

Description:  This script is a pipeline for quality-controlling transcriptome files and removing rRNA-like reads.  It uses mothur to eliminate any reads with homopolymers 8 bp or greater in length or any ambiguous bases.  It then uses mothur to dereplicate the data set, and then uses ribopicker to pull out rRNA sequences that match anythingi n the non-redundant rRNA database rrnadb (see the ribopicker website for how this database is assembled).  Finally, metaprocess.pl also uses ribopicker to select 16S/18S-like sequences from the removed rRNA sequences to facilitate taxonomic analysis based on small subunit rRNA.

Dependencies:  metaprocess.pl requires mothur and ribopicker as well as the rrnadb and ssr rRNA databases.  It also uses the POSIX, Getopt::Std, and File::Path perl modules.

Syntax:
	perl metaprocess.pl <transcriptome reads fasta file> <logfile> <number of processors to use in mothur operations> <optional mothur groups file>


blastseqs.postblast.pl

Description: This script takes a blast output file and finds the hits in the original, blasted sequence file and outputs them as a new fasta file. It also checks the mothur names file created by metaprocess.pl to see if any of the hits were dereplicated and if so, pulls out each sequence from the parent file (to make sure the count of transcript hits is accurate).  Finally, blastseqs.postblast.pl also translates the entire hit sequence in the same frame identified during blast analysis and output this as a separate amino acid fasta file.

Note that this file was extensively reworked during development and contains a lot of additional functionality currently commented out.  The current version is the version that was used in the final analyses in our published work but we kept the code because it might be useful for future analyses.

Dependencies: Commented-out sections require BLAST and MOTHUR.  This script uses the POSIX perl module.

	perl blastseqs.postblast.pl <BLAST csv output file> <fasta file that was BLASTed> <KEGGur taxonomy file> <cutoff % of sequence that must match query, 0-1> <output directory> <names file that matches BLASTed file>

UPDATE:
Some modified scripts were written to handle mothur data that uses a count table instead of names and groups files.  Also, the keggur.pl script was split into two sections -- the "pulldown" section which gathers sequences from KEGG, and the rest of the script, which aligns sequences and produces blast databases.  This was done to accomodate the insertion of user-supplied sequences to the KEGG sequences.

keggur.pulldown.pl
keggur.man.pl

These scripts take the same arguments as the original keggur.pl script.  Each is 1/2 of the keggur.pl script.  keggur.pulldown.pl stops after the KEGG database sequences are collected from the internet.  keggur.man.pl proceeds withe the alignment and other processing steps.  This allows introduction of user-curated sequences into the KEGG fasta file (and taxonomy files) prior to alignment.

blastseqs.count.postblast.pl

As blastseqs.postblast.pl, but takes a count_table instead of a groups file

Syntax:
	perl blastseqs.count.postblast.pl <BLAST csv output file> <fasta file that was BLASTed> <KEGGur taxonomy file> <cutoff % of sequence that must match query, 0-1> <output directory> <count_table that matches BLASTed file>

manalign.count.pl

Syntax exactly the same as manalign.pl; only change is the way that sequences are named in the output file, to facilitate use of the copyseqs.pl script to introduce copy numbers from the count table.

fixnames.mod.pl

Syntax as the fixnames.pl file. Modifications:
1. "tab" character appropriate for LINUX as opposed to UNIX used
2. Option to NOT use taxonomy file introduced (strange omission from original program)

copyseqs.pl

This program inserts replicate copies of a sequence into a fasta file based on counts within a count table file.  Note, this only works right now for one-group count table files and will blow up if you try to use more complex count tables. The purpose of this script is to make a mothur count table analyzed file amenable to quantitative analysis in pplacer. A superior script would directly modify the pplacer .jplace file, but I have not yet had time to figure out how to make that happen.

Syntax:
	perl copyseqs.pl <Sequence fasta file> <Output directory> <count table file that matches the sequence file>
