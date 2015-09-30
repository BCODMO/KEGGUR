#!/usr/local/bin/perl
# takes fasta sequence file and produces a mothur groups file
# divvying up sequences into groups
# arg0 = sequence filename
# arg1 = regexp indicating how to find group name
# arg2 = output filename
# arg3 = logfile

use strict;
use warnings;
use POSIX;
use Getopt::Std;
use File::Path qw( make_path );

# KEGGUR version number
my $version = 0.1;

my $date = POSIX::strftime("%A, %B %d, %Y",localtime());

#initialize variables and assign defaults if necessary

my ($filename,$regexp,$outputfile,$logfile)=@ARGV;
my ($inhandle,$outhandle);
my %lookup=();

unless ($logfile) {$logfile="makegroups.logfile"}
open (my $handle, ">>", $logfile);
print $handle "KEGGUR v $version\nmakegroups\n$date\n\n";
print $handle "Extracting mothur groups file from $filename\n";

unless (open($inhandle,$filename)) {
	print $handle "No such input file\n";
	print "No such input file\n";
	die
}

open ($outhandle,">$outputfile");

while(my $line = <$inhandle>) {
	unless ($line =~ /^>/) {next}
	chomp($line);
	my $group = $line;
	$group =~ s/$regexp/$1/;
	$line =~ s/>(\S*) .*$/$1/;
	print $outhandle "$line\t$group\n";
}

close ($inhandle);
close ($outhandle);

if (-e $outputfile) {
	print $handle "\nWrote new groups file, using the regular expression $regexp, to $outputfile\n\n";
}
else {
	print $handle "\nFailed to create groups file from $filename\nOperation aborted\n";
	print "\nFailed to create groups file from $filename\nOperation aborted\n";
	die
}

print $handle "KEGGUR makegroups complete\n";
close $handle;

exit;
