#!/usr/local/bin/perl
# Uses blastx to find matches to proteins in a KEGGur database in a sequence file
# takes the blast hits .csv file with the first value = hit sequence
# # and retrieves the full sequence record from the sequence file
# # Screens out hits that don't cover $cutoff % of parent read length
# # Also pull a names file for mothur
# # Last, gets taxonomy of best hit from master taxonomy file
# # Then classifies sequences using mothur wang classifier
# # Arg0 = BLAST CSV output file
# # Arg1 = Original blasted sequence file
# Arg 2 = Keggur taxonomy file
# Arg3 = cutoff (0-1) for alignment length to keep a hit
# Arg 4 = output directory
# Arg 5 = Original blasted sequence COUNT TABLE file

use strict;
use warnings;
use POSIX;

# KEGGUR version number
my $version = 0.1;

my $date = POSIX::strftime("%A, %B %d, %Y",localtime());

#initialize variables and assign defaults if necessary

#CHANGE OUTPUT NAMES TO OUTPUTDIR/SEQNAME.KEGGDB. etc etc
my ($csv,$seq,$taxfile,$cutoff,$outputdir,$countfile)=@ARGV;
my @csv =();
my (%lookup,%tax,%besthit,%frame)=();
my $flag=0;
my ($handle,$outhandle,$aaouthandle);
#my $seqname = $seq;
#my $db = $keggdb;
#$seqname =~ s/^.*\///;
#$db =~ s/^.*\///;



unless ($outputdir) {$outputdir = "./"}
unless (-d $outputdir) {system("mkdir $outputdir")}

my $filename=$csv;
$filename =~ s/^.*\///;
$filename="$outputdir/$filename";
my $logfile = "$filename.postblast.logfile";

# if (-f $logfile){
# 	print "Analysis has already been performed; cancelling.\n";
# 	die
# }

open (my $log, ">",$logfile);
print $log "KEGGUR v $version\nblastseqs.postblast\n$date\n\n";
print "KEGGUR v $version\nblastseqs.postblast\n$date\n\n";
print $log "Processing BLAST Output file $csv\n\n";
print "Processing BLAST Output file $csv\n\n";

unless (-e $taxfile) {
	print $log "Taxonomy file doesn't exist\nOperation aborted\n\n";
	print "Taxonomy file doesn't exist\nOperation aborted\n\n";
	die
}

unless (-e $countfile) {
	print $log "Count file doesn't exist\nOperation aborted\n\n";
	print "Count file doesn't exist\nOperation aborted\n\n";
	die
}

# unless ($groupnames) {$groupnames=$seqname}

# unless ($groupsfile) {
# 	system("mothur \"#set.logfile(name=$logfile, append=T);make.group(fasta=$seq,groups=$groupnames)\"");
# 	$groupsfile=$seq;
# 	$groupsfile =~ s/\.[^.]*$//;
# 	$groupsfile.=".groups";
# }

if ($cutoff > 1 or $cutoff < 0) {
	print $log "Cutoff value must be between 0 and 1.\n";
	print "Cutoff value must be between 0 and 1.\n";
	die
}

# if ($classifycutoff > 100 or $classifycutoff < 50) {
# 	print $log "Classification cutoff value must be between 50% and 100%.\n";
# 	print "Classification cutoff value must be between 50% and 100%.\n";
# 	die
# }
# 
# system("blastx -version 1>>$logfile 2>>$logfile");
# 
# system("blastx -best_hit_overhang=.25 -best_hit_score_edge=0.05 -num_threads $processors -db $keggdb.aa.db -query $seq -out $filename.csv -max_target_seqs 1 -outfmt \"10 qseqid sseqid qseq sseq evalue bitscore length gaps pident qcovs qstart qend sstart send qframe mismatch\" -evalue=0.0001 1>>$logfile 2>>$logfile");

unless (open($handle,"$csv")) {
	print $log "\nError opening BLAST output file\n";
	print "\nError opening BLAST output file\n";
	die
}


# if (-z "$filename.csv") {
# 	print $log "\nNo $keggdb BLAST hits in $seq\n";
# 	print $log "KEGGUR BLAST processing complete\n";
# 	print "\nNo $keggdb BLAST hits in $seq\n";
# 	print "KEGGUR BLAST processing complete\n";
# 	close $log;
# 	die
# }

@csv=<$handle>;
close($handle);

my $numrecords = scalar @csv;

# Make a taxhash

open($handle,$taxfile);
while (my $line = <$handle>) {
        chomp ($line);
        my $name = $line;
        my $lineage= $line;
        $name =~ s/\t.*$//;
        $lineage =~ s/^.*\t//;
        $tax{$name}=$lineage;
}
close($handle);

#process csv and combine multiple alignments from same read
# Removed this functionality because IT'S ONLY GETTING THINGS IN DIFFERENT FRAMES!
# for (my $i = 0; $i < ((scalar @csv)-1); $i++) {
#         my $name = $csv[$i];
#         $name =~ s/,.*$//;
#         my $nextname = $csv[$i+1];
#         $nextname =~ s/,.*$//;
#         if ($name eq $nextname) {
#                 my $sequence = $csv[$i];
#                 chomp($sequence);
#                 $sequence =~ s/^[^,]*,[^,]*,([^,]*),.*$/$1/;
#                 my $nextsequence = $csv[$i+1];
#                 chomp($nextsequence);
#                 $nextsequence =~ s/^[^,]*,[^,]*,([^,]*),.*$/$1/;
#                 $sequence .= $nextsequence;
#                 $csv[$i+1] =~ s/^([^,]*,[^,]*,)[^,]*/$1$sequence/;
#                 $csv[$i] =~ s/^([^,]*,[^,]*,)[^,]*/$1/;
#         }
# }

# process csv and make a hash of queries vs alignment lengths, queries vs. subjects, queries vs frames
foreach my $element (@csv) {
	chomp ($element);
	my $read = $element;
	$read =~ s/,.*$//;
	my $sequence = $element;
	$sequence =~ s/^[^,]*,[^,]*,([^,]*),.*$/$1/;
	# Eliminate gaps
	$sequence =~ s/-//g;
	# Eliminate anything after a stop codon
	# Note suspending this functionality for now; passing no judgments about pseudogenes
	#$sequence =~ s/\*.*$//;
	my $length = length($sequence);
	$lookup{$read}=$length;
	my $besthit = $element;
	$besthit =~ s/^[^,]*,([^,]*),.*$/$1/;
	unless (exists $tax{$besthit}) {
		print $log "$besthit not in taxonomy file\n\n";
		print "$besthit not in taxonomy file\n";
		die}
	$besthit{$read}=$besthit;
	# Make a hash of frames for each hit
	my $qframe = $element;
        $qframe =~ s/^[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,([^,]*),.*$/$1/;
        $frame{$read}=$qframe;
}

close($handle);

#get full sequences of any hits that clear the cutoff and make new nt and translated aa fasta files

my $seqhandle;

unless (open($seqhandle,$seq)) {print "No such seq file\n";die}

my $output = "$filename.nt.fasta";
my $aaoutput = $filename.".aa.fasta";
unless (open($outhandle,">$output")) {print "Can't open output file\n";die}
unless (open($aaouthandle,">$aaoutput")) {print "Can't open output file\n";die}

my $short = 0;
my $count = 0;
my @records = ();
while (my $record=<$seqhandle>) {
        chomp($record);
        unless ($record =~ />/) {next}
        my $name = $record;
        $name =~ s/^>(\S*).*$/$1/;
        if (exists $lookup{$name}) {
			if ($lookup{$name} eq "1") {print "Duplicate record $name\n";die}
			my $hitsequence = <$seqhandle>;
			chomp ($hitsequence);
			my $length = length ($hitsequence);
			my $limit = $length * $cutoff;
			my $lookupname = $lookup{$name} *3;
			unless ($lookupname >= $limit) {
				print $log "\n$name not long enough: $lookupname bp out of $length bp\n";
				delete $lookup{$name};
				delete $besthit{$name};
				$short++;
				next;
			}
			my $hitframe = $frame{$name};
			my $transsequence = translateseq($hitsequence,$hitframe);
			# Look for record in count file
			$count = `grep -P "$name\t" $countfile`;
			chomp $count;
			$count =~ s/^.*\t//;
			for (my $prk=1; $prk <= $count; $prk++) {
				print $outhandle ">$name"."_$prk\n";
				print $outhandle "$hitsequence\n";
				print $aaouthandle ">$name\n$transsequence"."_$prk\n";
				print $log "\nWrote $name to $output and $aaoutput\n";
			}
#			my $copies = scalar @records;
# 			if ($copies > 1) {
# 				print $log "\n$copies copies of $name found in $namesfile\n";
#                print "\n$copies copies of $name found in $namesfile\n";
# 			}
			$lookup{$name}="1"
        }
}

print $log "\n$short of $numrecords alignments were too short\n";
print "\n$short of $numrecords alignments were too short\n";


close ($outhandle);
close ($seqhandle);

# Get taxonomy of best BLAST hit and write to $filename.besthit.tax

my $hittax = $filename;
$hittax .= ".besthit.tax";
open ($handle,">$hittax");
while (my ($key,$value) = each(%besthit)) {
        my $whoosit = $besthit{$key};
        my $lineage = $tax{$whoosit};
        print $handle "$key\t$lineage\n";
        }
close($handle);

# Make mothur names files for blast hits

# if (open($handle,$namesfile)) {
#         my $namesfileout = "$filename.names";
#         open (my $namesout,">$namesfileout");
#         while (my $record=<$handle>) {
#                 chomp($record);
#                 my $name = $record;
#                 $name =~ s/^(\S*).*$/$1/;
#                 if (exists $lookup{$name}) {
#                         print $namesout "$record\n";
#                         my $morenames = $record;
#                         $morenames =~ s/^.*\t//;
#                         my @morenames = split(',',$morenames);
#                         foreach my $element (@morenames) {
#                                 $lookup{$element}=1;
#                         }
#                 }
#         }
#         close($handle);
#         close ($namesout);
# }elsif($namesfile){print "Can't open names file\n";die}

# if (open ($handle,$groupsfile)) {
#         my $groupsfileout="$filename.groups";
#         open (my $groupsout,">$groupsfileout");
#         while (my $record=<$handle>) {
#              chomp($record);
#              my $name = $record;
#              $name =~ s/^(\S*).*$/$1/;
#              if (exists $lookup{$name}) {
#                  print $groupsout "$record\n";
#              }
#         }
#         close($handle);
#         close ($groupsout);
# }elsif($groupsfile){print "Can't open groups file\n";die}

# Use mothur to classify sequences
# Note turning off this functionality for now; JUST processing blast hits file

#system("mothur \"#set.logfile(name=$logfile, append=T);classify.seqs(fasta=$filename.nt.fasta,name=$namesfile,cutoff=$classifycutoff,method=wang,template=$keggdb.template.fasta,taxonomy=$taxfile, processors=$processors,group=$groupsfile);summary.tax(taxonomy=$filename.besthit.tax,group=$groupsfile,reftaxonomy=$taxfile,name=$namesfile)\"");
 
print $log "KEGGUR BLAST processing complete\n";

close $log;

exit;

#-------------------------------------------------------------------------------------
sub codon2aa {
	# converts ARG frame into an aa and returns it
	# From book Beginning Perl For Bioinformatics

use strict;
use warnings;


     my($codon,$handle) = @_;
     if ( $codon =~ /GC./i)        { return 'A' }    # Alanine
     elsif ( $codon =~ /TG[TC]/i)     { return 'C' }    # Cysteine
     elsif ( $codon =~ /GA[TC]/i)     { return 'D' }    # Aspartic Acid
     elsif ( $codon =~ /GA[AG]/i)     { return 'E' }    # Glutamic Acid
     elsif ( $codon =~ /TT[TC]/i)     { return 'F' }    # Phenylalanine
     elsif ( $codon =~ /GG./i)        { return 'G' }    # Glycine
     elsif ( $codon =~ /CA[TC]/i)     { return 'H' }    # Histidine
     elsif ( $codon =~ /AT[TCA]/i)    { return 'I' }    # Isoleucine
     elsif ( $codon =~ /AA[AG]/i)     { return 'K' }    # Lysine
     elsif ( $codon =~ /TT[AG]|CT./i) { return 'L' }    # Leucine
     elsif ( $codon =~ /ATG/i)        { return 'M' }    # Methionine
     elsif ( $codon =~ /AA[TC]/i)     { return 'N' }    # Asparagine
     elsif ( $codon =~ /CC./i)        { return 'P' }    # Proline
     elsif ( $codon =~ /CA[AG]/i)     { return 'Q' }    # Glutamine
     elsif ( $codon =~ /CG.|AG[AG]/i) { return 'R' }    # Arginine
     elsif ( $codon =~ /TC.|AG[TC]/i) { return 'S' }    # Serine
     elsif ( $codon =~ /AC./i)        { return 'T' }    # Threonine
     elsif ( $codon =~ /GT./i)        { return 'V' }    # Valine
     elsif ( $codon =~ /TGG/i)        { return 'W' }    # Tryptophan
     elsif ( $codon =~ /TA[TC]/i)     { return 'Y' }    # Tyrosine
     elsif ( $codon =~ /TA[AG]|TGA/i) { return '*' }    # Stop
     else {
         #print $handle "Bad codon \"$codon\"!!\n";
         return "X";
     }
 }

#-------------------------------------------------------------------------------------
sub translateseq {
# converts ARG0 fasta DNA sequence into peptide in ARG1 reading frame

use strict;
use warnings;


my ($seq,$frame,$handle)=@_;

my $pep='';
chomp($seq);

# If frame is -1, -2, or -3, get reverse compliment of sequence

if ($frame < 0) {
	$seq = revcom($seq);
	$frame = abs($frame);
}

# If frame is -2, -3, 2, or 3, remove 1 or 2 bp from first to put in proper frame.

$seq = substr($seq,$frame-1);

# Remove bp from the end if they aren't a multiple of 3

my $inframe = length($seq);
if (my $mod=$inframe % 3) {
	$seq=substr($seq,0,$inframe-$mod);
}

# Do actual translation

for (my $ntp=0; $ntp < (length($seq)-2); $ntp=$ntp+3) {
			my $codon=substr($seq,$ntp,3);
			my $aa=codon2aa($codon,$handle);
		# If stop codon, stop translating (eliminating this for now)
# 		if ($aa eq "*") {
# 			$pep=$pep.$aa;
# 			last
# 		}
#		else {
			$pep=$pep.$aa;
#		}
	}

return $pep;
}

#--------------------------------------------------------------------------------
sub revcom {

# Code from Beginning Perl For Bioinformatics

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}

#---------------------------------------------------------------------------------
