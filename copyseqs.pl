#!/usr/local/bin/perl
# Replicates a sequence in a fasta file based on copy number in a count table
# # Arg0 = Sequence file
# Arg 1 = output directory
# Arg 2 = COUNT TABLE file

use strict;
use warnings;
use POSIX;

# KEGGUR version number
my $version = 0.1;

my $date = POSIX::strftime("%A, %B %d, %Y",localtime());

#initialize variables and assign defaults if necessary

#CHANGE OUTPUT NAMES TO OUTPUTDIR/SEQNAME.KEGGDB. etc etc
my ($seq,$outputdir,$countfile)=@ARGV;
my @csv =();
my (%lookup,%tax,%besthit,%frame)=();
my $flag=0;
my ($handle,$outhandle,$aaouthandle);

unless ($outputdir) {$outputdir = "./"}
unless (-d $outputdir) {system("mkdir $outputdir")}

my $filename=$seq;
$filename =~ s/^.*\///;
$filename="$outputdir/$filename";
my $logfile = "$filename.copyseqs.logfile";

# if (-f $logfile){
# 	print "Analysis has already been performed; cancelling.\n";
# 	die
# }

open (my $log, ">",$logfile);
print $log "KEGGUR v $version\ncopyseqs\n$date\n\n";
print "KEGGUR v $version\ncopyseqs\n$date\n\n";
print $log "Replicating sequences based on $countfile\n\n";
print "Replicating sequences based on $countfile\n\n";

unless (-e $countfile) {
	print $log "Count file doesn't exist\nOperation aborted\n\n";
	print "Count file doesn't exist\nOperation aborted\n\n";
	die
}

my $seqhandle;

unless (open($seqhandle,$seq)) {print "No such seq file\n";die}

my $output = "$filename.rep.fasta";
unless (open($outhandle,">$output")) {print "Can't open output file\n";die}

my $short = 0;
my @count = ();
my @records = ();
while (my $record=<$seqhandle>) {
        chomp($record);
        unless ($record =~ />/) {next}
        my $name = $record;
        $name =~ s/^>(\S*).*$/$1/;
		my $hitsequence = <$seqhandle>;
		chomp ($hitsequence);
			# Look for record in count file
			# my @number = grepcount($countfile,$name);
# 			foreach my $element (@number) {print chomp $element}
# 			if ($number[1]=0) {$count=1;print "yes\n"}
# # 			else {
				@count = `grep "$name" $countfile`;
				my $count = $count[0];
				if (!defined($count)) {$count = 1}
				$count =~ s/^.*\t//;
# 			}
			for (my $prk=1; $prk <= $count; $prk++) {
				if ($count > 1) {
					print $outhandle ">$name"."_$prk\n";
				}
				else {
					print $outhandle ">$name\n";
				}
				print $outhandle "$hitsequence\n";
				print $log "\nWrote $name to $output \n";
			}
}

close ($outhandle);
close ($seqhandle);
 
print $log "KEGGUR copy replication complete\n";

close $log;

exit;

#-----------------------------------------------

sub grepcount{
#takes 2 args, 0 = file name, 1 = string to search for
#uses UNIX grep -c to count # of lines containing string
#returns an array of numbers

use strict;
use warnings;

my ($name,$string) = @_;
my @count=();

my @test =`grep -c \"$string\" $name`;
foreach my $element (@test) {
        chomp $element;
        $element =~ s/^[^:]*://;
        $element .= "\n";
        push (@count,$element);
}

return @count;
}
