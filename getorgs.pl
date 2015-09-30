
# produces a master taxonomy of all the organism refesqs in the kegg database

use strict;
use warnings;
use POSIX;

# KEGGUR version number
my $version = 0.1;

my $date = POSIX::strftime("%A, %B %d, %Y",localtime());

open (my $handle, ">","getorgs.logfile");
	print $handle "KEGGUR v $version\ngetorgs\n$date\n\n";
close $handle;

getorgs();
buildtaxonomy('kegg.orgs','kegg.tax');

exit;

#------------------------------------------------------------------------
sub getorgs {
#downloads a complete list of codes for refseqs in the KEGG database
#takes output filename as an arg or defaults to 'kegg.orgs'

use strict;
use warnings;


my ($outputfile)=@_;

unless ($outputfile) {$outputfile='kegg.orgs'}

my @list=();

system ('wget -O '.$outputfile.'.temp http://www.genome.jp/kegg-bin/show_organism?category=* -a getorgs.logfile');

open(my $handle,"$outputfile.temp");
	@list=<$handle>;
close ($handle);

system("rm $outputfile.temp 1>>getorgs.logfile 2>>getorgs.logfile");

@list=grep(/show_organism\?org=/,@list);

foreach my $element (@list) {
	$element =~ s/^.*show_organism[^>]*>//g;
	$element =~ s/<.*$//g;
}

unless (open($handle,'>'.$outputfile)) {die "Can't open file1\n"}
	foreach my $element (@list) {
		print $handle $element;
	}
close($handle);

open ($handle,">>", "getorgs.logfile");
	print $handle "\nRetrieved current organism file list from KEGG database\n";
	print $handle "Saved file list as $outputfile\n\n";
close $handle;

return @list;
}

#----------------------------------------------------------

sub buildtaxonomy {
# takes list of KEGG organism codes and builds a mothur-compatible
# taxonomy file. Takes input and output filenames as args 0 and 1.
#
use strict;
use warnings;


my ($inputfile,$outputfile,$namesfile)=@_;
my @list=();

unless ($outputfile) {$outputfile="$inputfile.tax"}

if(open (my $handle, $inputfile)) {
        @list=<$handle>;
        close($handle);
}
else {
        die "Can't open file $inputfile\n";
        exit;
}

open (my $header,'>'.$outputfile);

foreach my $element (@list) {
	my ($name,$taxonomy)='';
	chomp $element;
        system('wget -O '.$inputfile.'.temp http://www.genome.jp/dbget-bin/www_bget?'.$element.' -a getorgs.logfile');
        open(my $handle,"$inputfile.temp");
                my @stuff=<$handle>;
        close($handle);
        system("rm $inputfile.temp 1>>getorgs.logfile 2>>getorgs.logfile");
        my $flag=0;
	SCAN:
        foreach my $line (@stuff) {
                if ($line =~ ">Definition<") {
                        $flag=1;
                }
                elsif ($line =~ "&nbsp;&nbsp;Lineage") {
                        $flag=2;
                }
                else {
                        if ($flag==1) {
                                $line =~ s/<[^>]*>//g;
				$line =~ s/;//g;
                                chomp $line;
                                $name=$line;
                                $flag=0;
                                next;
                        }
                        elsif ($flag==2) {
                                $line =~ s/<[^>]*>//g;
                                $line =~ s/;\s*/;/g;
                                chomp $line;
				$line =~ s/\s/_/g;
                                $taxonomy=$line.';';
                                $flag=0;
                                last SCAN;
                        }
                        else {next};
                }
        }
	unless ($name && $taxonomy) {print $header "Insufficient data for $element\n";die}
	my $newline = $element.'|'.$name."\t$taxonomy\n";
	print $header $newline;
}

close $header;

open (my $handle, ">>","getorgs.logfile");
	print $handle "\nRetrieved taxonomic information for all current KEGG organisms\n";
	print $handle "Saved as $outputfile\n\n";
close $handle;

return;

}
