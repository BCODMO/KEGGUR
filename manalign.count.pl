#!/usr/local/bin/perl
#Manual Alignment based on Blast results vs. master alignment
#Arg0 = database reference file (nt's)
#Arg1 = postblast hit nt file
#Arg2 = BLAST csv file (need sseq, qseq, sstart, qstart, qframe)
#Arg3 = output file
# Modified for use with a count table; very small change in how sequences are named when written to output file

use strict;
use warnings;
use POSIX;

# KEGGUR version number
my $version = 0.1;

my $date = POSIX::strftime("%A, %B %d, %Y",localtime());

#initialize variables

my (@csv,@qalign,@qcodes) = ();
my (%ref,%hits,%sseqid,%qseq,%sseq,%qstart,%sstart,%qframe,%qalign,%qcodes)=();

my ($ref,$hits,$csv,$output)=@ARGV;

my $handle;

my $n=0;
%ref=readref($ref);

#Make sure all aligned reference sequences are the same length and get that length
my $reflength=0;
foreach (sort keys %ref) {
	my $thislength = length $ref{$_};
	if ($reflength==0) {
		$reflength=$thislength;
		next
	}
	if ($reflength != $thislength) {
		print "Reference alignment sequences are not all the same length\n";
		die
	}
}


%hits=readfasta($hits);

unless (open($handle,"$csv")) {
 #	print $log "\nError opening BLAST output file\n";
 	print "\nError opening BLAST output file\n";
	die
}
@csv=<$handle>;
close($handle);

#prepare hashes from blast csv

foreach my $element (@csv) {
	chomp ($element);
	
	my $qseqid=$element;
	my $sseqid=$element;
	my $qseq=$element;
	my $sseq=$element;
	my $qstart=$element;
	my $sstart=$element;
	my $qframe=$element;
	$qseqid =~ s/,.*$//;	
	$sseqid=~s/^[^,]*,([^,]*),.*$/$1/;
	$qseq=~s/^[^,]*,[^,]*,([^,]*),.*$/$1/;
	$sseq=~s/^[^,]*,[^,]*,[^,]*,([^,]*),.*$/$1/;	$sstart=~s/^[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,([^,]*),.*$/$1/;	$qframe=~s/^[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,([^,]*),.*$/$1/;
	if ($qframe < 0) {
		$qstart=~s/^[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,([^,]*),.*$/$1/;
	}
	else {
		$qstart=~s/^[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,[^,]*,([^,]*),.*$/$1/;
	}
	$sseqid{$qseqid}=$sseqid;
	$qseq{$qseqid}=$qseq;
	$sseq{$qseqid}=$sseq;
	$qstart{$qseqid}=$qstart;
	$sstart{$qseqid}=$sstart;
	$qframe{$qseqid}=$qframe;
}

#produce initial aligned sequences; get length of longest one

my $longestsequence = 0;

foreach my $qseqid (@csv) {
	# check if the qseq is in the hits fasta, if not move on
	$qseqid =~ s/,.*$//;
	if (exists $hits{$qseqid}) {
		#get hitseq, check frame and reverse if necessary
		my $hitseq=$hits{$qseqid};
		my $checkframe=$qframe{$qseqid};
		if ($checkframe < 0) {$hitseq=revcom($hitseq)}
		#get alignment and code
		my $xseq = "$sseqid{$qseqid}";
		my @align = manalign($ref{$sseqid{$qseqid}},$hitseq,$sseq{$qseqid},$qseq{$qseqid},$sstart{$qseqid},$qstart{$qseqid});
		if (length($align[1]) > $longestsequence) {
			$longestsequence = length($align[1]);
		}
		$qalign{$qseqid}=$align[0];
		$qcodes{$qseqid}=$align[1];
	}
}

# Merge refalignment w. aligned hits
@qalign{keys %ref} = values %ref;

#Add refs into qcodes as all dashes, no pluses

  foreach (sort keys %ref) {
     $qcodes{$_}="-"x$reflength;
   }

# go through and insert missing gaps

for (my $search = 0; $search < $longestsequence; $search++) {
	my $plus="";
	foreach (sort keys %qcodes) {
		if (length $qcodes{$_} < $search) {next}
		$plus .= substr($qcodes{$_},$search,1);
	}
	if (index($plus,"+") != -1) {
		#There's at least one new refseq gap required at this site; add gap to everything without a + in this position in its code
		foreach (sort keys %qalign) {
 			if (substr($qcodes{$_},$search,1) eq "-") {
 				substr($qalign{$_},$search,0)="-";
			}
		}
	}
}

# Check to make sure all are the same length
my $newlength=0;
foreach (sort keys %qalign) {
	my $thislength = length $qalign{$_};
	if ($newlength==0) {
		$newlength=$thislength;
		next
	}
	if ($newlength != $thislength) {
		print "New alignment sequences are not all the same length\n";
		print "$_\n";
		die
	}
}

print "Original alignment was $reflength bp long; new alignment is $newlength bp long.\n";

# Write new alignment file

my $alignref = \%qalign;
writefasta($alignref,$output);

exit;

#-------------------------------------
sub manalign {

#"manual" alignment
#arg0=refseq (in nt's, from main alignment)
#arg1=hitseq (in nt's, from main metagenome)
#arg2=sseq (in aa's, from blast csv)
#arg3=qseq (in aa's, from blast csv)
#arg4=sstart
#arg5=qstart

#initialize variables and assign defaults if necessary

my ($refseq,$hitseq,$sseq,$qseq,$sstart,$qstart)=@_;

my (@alignment) = ();

my $alength = length($sseq);

if ($alength != length($qseq)) {
	print "Alignments aren't the same length; aborting\n";
	die
}

#in real iteration, will take seqids and find original nt files
#still needs to incorporate reversing negative frames, right now only works when qend>qstart

my ($outseq,$codeseq,$refchar,$hitchar,$saa,$qaa) = "";
my $i = 0;
my $qcount=0;
my $acount=0;

my $qshift=($sstart-1)*3-$qstart+1;
my $reflength = length($refseq);
if ($reflength < length($hitseq)) {
	print "Target sequence $hitseq is longer than the reference sequence $refseq; skipping\n";
	print "Sseq = $sseq\nQseq = $qseq\n";
}

for ($i=0; $i < $reflength; $i++) {
	$refchar = substr($refseq,$i,1);
	#skip blanks in reference alignment
	if ($refchar eq "-") {
		$outseq.="-";
		$codeseq.="-";
		next
	}
	#add leading blanks up to beginning of where sequence starts
	if ($qcount < $qshift) {
		$outseq.="-";
		$codeseq.="-";
		$qcount++;
		next
	}
	#until alignment start, add lower-case letters from hitseq aligned exactly with refseq
	if ($qcount < ($sstart-1)*3) {
		$outseq.=lc(substr($hitseq,$qcount-$qshift,1));
		$codeseq.="-";
		$qcount++;
		next
	}

# 	#during alignment, check each alignment position and add in codons (groups of 3), uppercase
	if ($qcount < ($sstart-1)*3+length($sseq)*3) {
		$saa=substr($sseq,$acount,1);
		$qaa=substr($qseq,$acount,1);
		if ($qaa eq "-") {
			$outseq.="-"x3;
			$codeseq.="-"x3;
			#implement something to make sure ref doesn't have blanks in middle of codons
			$acount++;
			$qcount+=3;
			$qshift+=3;
			$i+=2;
			next
		}
		if ($saa eq "-") {
			#add +'s if gaps in the ref alignment and don't increment the main counter
			$codeseq.="+"x3;
			$outseq.=substr($hitseq,$qcount-$qshift,3);
			$acount++;
			$qcount+=3;
			$i-=1;
			next
		}
		$outseq.=substr($hitseq,$qcount-$qshift,3);
		$codeseq.="-"x3;
		$qcount+=3;
		$acount++;
		$i+=2;
		next
	}
	
	# After alignment add remaining sequence (lower case) and pad to end of refseq with blanks
	if ($qcount < $qshift+length($hitseq)) {
		$outseq.=lc(substr($hitseq,$qcount-$qshift,1));
		$codeseq.="-";
		$qcount++;
		next
	}
	$outseq.="-";
	$codeseq.="-";
}

# Go from front to back of aligned sequences. If any have a "+" in codeseq, insert gap in refseq
# for ($i=0; $i < length($codeseq); $i++) {
# 	if (substr($codeseq,$i,1) eq "+") {substr($refseq,$i,0)="-"}
# }

 # print "$refseq\n";
 # print "$outseq\n";
#  
 #print "$codeseq\n";
 #die;

@alignment = ($outseq,$codeseq);

return @alignment;

}

#----------------------------------------------------------------------------------------------------------
sub readfasta # return hash of headers and sequences from file ARG0. modified for this program to kill everything but name
{

use strict;
use warnings;


my ($file)=@_;

my $sequence='';
my $line='';
my $name='';
my %fasta = ();

open(my $handle,$file);
	my @data=<$handle>;
close($handle);

foreach my $line (@data) {

       if ($line =~ /^\s*$/) {
            next;

       } elsif($line =~ /^\s*#/) {
            next;

       } elsif($line =~ /^>/) {
	    $sequence =~ s/\s//g;
	    if ($name) {
	    	chomp $sequence;
	    	$fasta{$name}=$sequence;
		$sequence='';
	    }
	    $name=$line;
	    $name =~ s/>//;
	    $name =~ s/\s.*$//;
	    $name =~ s/_[0-9]*$//;
	    chomp $name;
	
       } else {
	    chomp($line);
            $sequence .= $line;
       }

$fasta{$name}=$sequence;

}


return %fasta;
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

#------------------------------------------------------------------------------------

sub writefasta # takes referenced hash of headers and sequences ARG0 and writes a file or returns an array
{

use strict;
use warnings;


my ($fasta,$outfile)=@_;
my @fasta = ();

while (my ($name,$seq) = each (%{$fasta})) {
	chomp $name;
	chomp $seq;
	$name = ">".$name."\n";
	$seq = "$seq\n";
	push(@fasta,$name);
	push(@fasta,$seq);
}

if ($outfile) {
	open (my $handle,'>'.$outfile);
		foreach my $element (@fasta) {
			print $handle "$element";
		}
	close ($handle);
}

return @fasta;
}
#----------------------------------------------------------------------------------------------------------
sub readref # return hash of headers and sequences from file ARG0. modified for this program to kill everything but name
{

use strict;
use warnings;


my ($file)=@_;

my $sequence='';
my $line='';
my $name='';
my %fasta = ();

open(my $handle,$file);
	my @data=<$handle>;
close($handle);

foreach my $line (@data) {

       if ($line =~ /^\s*$/) {
            next;

       } elsif($line =~ /^\s*#/) {
            next;

       } elsif($line =~ /^>/) {
	    $sequence =~ s/\s//g;
	    if ($name) {
	    	chomp $sequence;
	    	$fasta{$name}=$sequence;
		$sequence='';
	    }
	    $name=$line;
	    $name =~ s/>//;
	    $name =~ s/\s.*$//;
	    chomp $name;
	
       } else {
	    chomp($line);
            $sequence .= $line;
       }

$fasta{$name}=$sequence;

}


return %fasta;
}
