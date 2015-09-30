#!/usr/local/bin/perl
# Pre-process a metagenome/metatranscriptome file for use with KEGGUR and mothur
# Requires mothur and ribopicker to be in the path
# Arg 0 = readsfile
# Arg 1 = logfile
# Arg 2 = processors
# Arg 3 = groupsfile
# If groups are in separate files and you want to analyze together, use mother make.group command.

use strict;
use warnings;
use POSIX;
use Getopt::Std;
use File::Path qw( make_path );

# KEGGUR version number
my $version = 0.1;

my $date = POSIX::strftime("%A, %B %d, %Y",localtime());

#initialize variables and assign defaults if necessary

my ($readsfasta,$logfile,$processors,$groups)=@ARGV;
my $newgroups;
unless ($logfile) {$logfile="metaprocess.logfile"}

if (!-e $readsfasta) {die "$readsfasta does not exist\n"}

if ($groups) {
	if (!-e $groups) {die "$groups does not exist\n"}
}

my $file = $readsfasta;
$file =~ s/.fa[^.]*$//;

open (my $handle, ">>", $logfile);
print $handle "KEGGUR v $version\nmetaprocess\n$date\n\n";
print $handle "Pre-processing $readsfasta for use with mothur/KEGGUR\n";

# Use mothur to eliminate bad reads

unless ($processors) {$processors=1}

my $goodfile = $readsfasta;
$goodfile =~ s/\.([^.]*)$/\.good\.$1/;
my $badfile = $readsfasta;
$badfile =~ s/\.[^.]*$/\.bad\.accnos/;
my $uniquefile = $goodfile;
$uniquefile =~ s/\.([^.]*)$/\.unique\.$1/;
my $names = $goodfile;
$names =~ s/\.[^.]*$/\.names/;

if ($groups) {
	system("mothur \"#set.logfile(name=$logfile, append=T);screen.seqs(fasta=$readsfasta,maxambig=0,maxhomop=8,group=$groups,processors=$processors);unique.seqs(fasta=$goodfile)\"");
	$newgroups = $groups;
	$newgroups =~ s/\.[^.]*$/\.good\.groups/;
}
else {
	$groups="$file.groups";
	system("mothur \"#set.logfile(name=$logfile, append=T);make.group(fasta=$readsfasta,groups=$file);screen.seqs(fasta=$readsfasta,group=$groups,maxambig=0,maxhomop=8,processors=$processors);unique.seqs(fasta=$goodfile)\"");
	$groups=$file;
	$newgroups="$groups.good.groups";
}

if (-e $goodfile) {
	my $badcount = `grep -c \"\" $badfile`;
	chomp $badcount;
	print $handle "\nUsed mothur to eliminate $badcount bad reads from $readsfasta.\n";
}
else {
	print "mothur screen.seqs failed\nOperation aborted\n";
	print $handle "\nmothur screen.seqs failed\nOperation aborted\n";
	die
}

if (!-e $uniquefile) {
	print $handle "\nmothur unique.seqs failed\nOperation aborted\n";
	print "\nmothur unique.seqs failed\nOperation aborted\n";
	die
}

my $preunique = `grep -c \">\" $goodfile`;
chomp $preunique;
my $postunique = `grep -c \">\" $uniquefile`;
chomp $postunique;
my $uniquecount = $preunique - $postunique;

print $handle "\nUsed mothur to eliminate $uniquecount duplicate sequences\n\n";

my $ribofile = $uniquefile;
my $ribodir= $uniquefile;
if (index($ribodir,"/") == -1) {
        $ribodir = './';
} else {
	$ribodir =~ s/\/[^\/]*$//;
	$ribofile =~ s/.*\///;
	unless ($ribodir) {$ribodir = './'}
}

$ribofile =~ s/.fa[^.]*$//;

# Use ribopicker to isolate non-rRNA reads

my $riboversion = `ribopicker.pl -version`;
print $handle "\n$riboversion\n\n";

system("ribopicker.pl -c 50 -i 90 -f $uniquefile -db rrnadb -id $ribofile -z 3 -out_dir $ribodir 1>>$logfile 2>>$logfile");

if (-e "$ribodir/".$ribofile."_nonrrna.fa") {
	my $ribocount = `grep -c \">\" $ribodir/$ribofile\_rrna.fa`;
	chomp $ribocount;
	print $handle "\nUsed riboPicker to identify $ribocount sequences with matches in the non-redundant rRNA database rrnadb\n";
} else {
	print $handle "\nriboPicker failed.\nOperation aborted\n";
	print "\nriboPicker failed.\nOperation aborted\n";
	die;
}

# Use ribopicker to isolate small subunit rRNA reads

system("ribopicker.pl -c 50 -i 90 -f $ribodir/".$ribofile."_rrna.fa -db ssr -id $ribofile.16S.silva -z 3 -out_dir $ribodir 1>>$logfile 2>>$logfile");

system("rm $ribodir/$ribofile.16S.silva_nonrrna.fa 1>>$logfile 2>>$logfile");

if (-e "$ribodir/$ribofile.16S.silva_rrna.fa") {
	my $ssucount = `grep -c \">\" $ribodir/$ribofile.16S.silva\_rrna.fa`;
	chomp $ssucount;
	print $handle "\nUsed riboPicker to identify $ssucount SSU rRNA sequences using SILVA database\n";
} else {
	print $handle "\nriboPicker failed.\nOperation aborted\n";
	print "\nriboPicker failed.\nOperation aborted\n";
	die;
}

# Remove rRNA reads and rename final files

system("mothur \"#set.logfile(name=$logfile, append=T);list.seqs(fasta=$ribodir/$ribofile\_rrna.fa);remove.seqs(accnos=$ribodir/$ribofile\_rrna.accnos,group=$newgroups);remove.seqs(accnos=$ribodir/$ribofile\_rrna.accnos,name=$names)\"");
my $pickgroups = $newgroups;
$pickgroups =~ s/\.([^.]*)$/\.pick\.$1/;
unless (-e $pickgroups) {
	print $handle "\nFailed to remove rRNA from groups file\nOperation aborted\n\n";
	print "\nFailed to remove rRNA from groups file\nOperation aborted\n\n";
	die
}
system("mv $pickgroups $file.processed.nonrrna.groups 1>>$logfile 2>>$logfile");

my $picknames = $names;
$picknames =~ s/\.([^.]*)$/\.pick\.$1/;

unless (-e "$picknames") {
	print $handle "\nFailed to remove rRNA from names file\nOperation aborted\n\n";
	print "\nFailed to remove rRNA from names file\nOperation aborted\n\n";
	die
}

system("mv $ribodir/$ribofile\_nonrrna.fa $file.processed.nonrrna.fasta 1>>$logfile 2>>$logfile");
system("mv $picknames $file.processed.nonrrna.names 1>>$logfile 2>>$logfile");

print $handle "\nFinal processed mothur-ready files written: $file.processed.nonrrna.fasta, $file.processed.nonrrna.groups, and $file.processed.nonrrna.names\n\n";

print $handle "KEGGUR metagenomic processing complete\n";
close $handle;

exit;
