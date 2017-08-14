# Replaces all kegg organism codes in a document with full names.
# $file = file to work on
# $outputfile = file to write to
# $namefile = kegg organism codes + names file (kegg.names)
# $delimitera and b = what's immediately before and after the organism code
# $taxonomy = if given, will use named taxonomy file to replace name with taxonomy followed by species

use strict;
use warnings;

my ($file,$outputfile,$namefile,$delimitera,$delimiterb,$taxonomy)=@ARGV;

my (%tax) = ();

unless (-e $file) {
	print "Can't open $file\n";
	die
}

unless (-e $namefile) {
	print "Can't open $namefile\n";
}

open (my $if,$file);
	my @file=<$if>;
close $if;

open(my $nf,$namefile);
	my @list=<$nf>;
close($nf);

if ($taxonomy) {
	%tax=taxhash($taxonomy);


foreach my $element (@list) {
	chomp $element;
	my $code = $element;
	my $name = $element;
	$code =~ s/	.*$//;
	$name =~ s/^.*	//;
	my $thistax = $tax{$code};
	#print "$code\t$name\t$thistax\n";
	my $left = '\\'.$delimitera.$code.$delimiterb;
	my $right = $delimitera.$thistax.$name.$delimiterb;
	foreach my $line (@file) {
		$line =~ s/$left/$right/g;
	}
}
}
else {
foreach my $element (@list) {
	chomp $element;
	my $code = $element;
	my $name = $element;
	$code =~ s/	.*$//;
	$name =~ s/^.*	//;
	#print "$code\t$name\t$thistax\n";
	my $left = '\\'.$delimitera.$code.$delimiterb;
	my $right = $delimitera.$name.$delimiterb;
	foreach my $line (@file) {
		$line =~ s/$left/$right/g;
	}
}
}

open (my $of,">$outputfile");
	foreach my $line (@file) {
		print $of $line;
	}
close ($of);

exit;

#-------------------------------------------

sub taxhash{
# opens a mothur taxonomy file OR array reference and makes a hash w. codes as keys
# # and taxonomy as values.  Takes inputfile OR array variable 
# as arg and returns a hash
# note: returns ONLY code; trims name etc separated by |

use strict;
use warnings;


my ($input) = @_;
my (@tax) = ();
my (%taxhash) = ();

if (-e $input) {
	#my $handle='';
	open (my $handle,$input);
	@tax = <$handle>;
	close($handle);
} else {@tax = @{$input}}

foreach my $line (@tax) {
       my $code = $line;
       my $lineage = $line;
       $code =~ s/\|.*$//g;
       chomp $code;
       chomp $lineage;
       $lineage =~ s/^.*\t//g;
       $lineage =~ s/[\; ]/_/g;
       $lineage =~ s/[^a-zA-Z0-9_]//g;
       $taxhash{$code}=$lineage;
}
return %taxhash;
}
