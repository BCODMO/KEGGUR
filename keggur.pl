#!/usr/local/bin/perl
# Requires MOTHUR, MUSCLE, AND BLAST+
# Given a kegg enzyme code and a filename, this batch file will pull down
# all annotated genes from the kegg database with that code, translate them,
# align them in MUSCLE until 3 subsequent runs are identical, and back-translate
# the amino acid alignment to the original nucleotides (with gaps).
# requires that the MUSCLE executable, mothur, and blast+ 
# be on the same path as this file.
# $enzyme is the kegg enzyme code (orthology group)
# $file is the base name to use in naming all the files this program generates
# $outputdir is the directory to save everything
# $muscleiterations is max number of muscle runs to try before giving up on convergence or finding a loop
# $taxonomy is the master taxonomy file, obtained by running getorgs.pl
# NOTE: will blow up if the getorgs.pl taxonomy is old; i.e. if there are gene hits that aren't represented in the reference taxonomy.
# $cluster variable tells whether or not to cluster based on cutoff value
# $cutoff is the cutoff value for OTU clustering; irrelevant unless $cluster is 1

use strict;
use warnings;
use POSIX;
use Getopt::Std;
use File::Path qw( make_path );

# KEGGUR version number
my $version = 0.1;

my $date = POSIX::strftime("%A, %B %d, %Y",localtime());

#initialize variables and assign defaults if necessary

my ($enzyme,$file,$outputdir,$muscleiterations,$taxonomy,$cluster,$cutoff)=@ARGV;

unless ($enzyme) {die "Need an orthology code K#####\nSyntax: keggur.pl ",'$enzyme $file $outputdir $muscleiterations $cluster $cutoff',"\n"}
if ($enzyme !~ /K[0-9][0-9][0-9][0-9][0-9]$/) {die "Need an orthology code K#####\nSyntax: keggur.pl ",'$enzyme $file $outputdir $muscleiterations $cluster $cutoff',"\n"; }
unless ($file) {$file = $enzyme}
unless ($outputdir) {$outputdir='./'}
if (!-d $outputdir) {
	make_path $outputdir or die "Output directory could not be made\n"
}
unless ($cluster) {$cluster=0}
unless ($cutoff) {$cutoff=1}
unless ($muscleiterations) {$muscleiterations=10}
unless ($taxonomy) {$taxonomy="kegg.tax"}
$file = $outputdir.'/'.$file;
if ($cluster !~ /^[01]$/) {die "Cluster flag must be 0 or 1.\n"}
if ($cutoff > 1 or $cutoff < 0) {die "Cutoff value must be between 0 and 1.\n"}
if ($muscleiterations <= 3) {die "Invalid number of muscle iterations; must be >= 3.\n"}

my $logfile = $file.'.keggur.logfile';
open (my $handle, ">>",$logfile);
print $handle "KEGGUR v $version\nkeggur\n$date\n\n";
print $handle "Preparing mothur and BLAST database files for KEGG enzyme $enzyme\n";
print $handle "Saving files to $outputdir\n";
print $handle "Using at most $muscleiterations alignment passes using MUSCLE\n";
if ($cluster==1) {print $handle "OTU clustering is on with OTU cutoff of $cutoff\n\n"}
else {print $handle "OTU clustering is off\n\n"}

# get list of all genes with the target orthology code from the kegg database and saves to $file.list

my $listcount = getlist($enzyme,"$file.list",$logfile);
if (-e "$file.list"){
	print $handle "\nRetrieved current list of $listcount KEGG genes in orthology group $enzyme and saved to $file.list\n\n";
} else {
	print $handle "\nFailed to retrieve list of KEGG genes\nOperation aborted\n\n";
	print "\nFailed to retrieve list of KEGG genes\nOperation aborted\n\n";
	die
}

# gets nucleotide sequences of all genes and saves them to $file.nt.fasta

my $keggcount=getkegg("$file.list","$file.nt.fasta",$logfile);

if (-e "$file.nt.fasta") {
	if ($keggcount == $listcount) {
		print $handle "\nRetrieved nucleotide sequences of all $listcount $enzyme genes and saved to $file.nt.fasta\n\n";
	} else {
		print $handle "\nExpected $listcount sequences by only receives $keggcount sequences\nOperation aborted\n\n";
		print "\nExpected $listcount sequences by only receives $keggcount sequences\nOperation aborted\n\n";
		die
	}
} else {
	print $handle "\nFailed to retrieve nucleotide sequences\nOperation aborted\n\n";
	print "\nFailed to retrieve nucleotide sequences\nOperation aborted\n\n";
	die
}

# removes redundant and ambiguous sequences, and creates a names file.

system('mothur "#set.logfile(name='.$logfile.',append=T);trim.seqs(fasta='.$file.'.nt.fasta,maxambig=0);unique.seqs(fasta='.$file.'.nt.trim.fasta)"');

if (-e "$file.nt.trim.fasta") {
	my $scrapcount = grepcount("$file.nt.scrap.fasta",">");
	print $handle "\nRemoved $scrapcount sequences with ambiguous bases\n\n";
} else {
	print $handle "\nTrim.seqs failed\nOperation aborted\n\n";
	print "\nTrim.seqs failed\nOperation aborted\n\n";
	die
}

if (-e "$file.nt.trim.unique.fasta") {
	my $trimcount = grepcount("$file.nt.trim.fasta",">");
	my $unicount = grepcount("$file.nt.trim.unique.fasta",">");
	my $removed = $trimcount - $unicount;
	print $handle "\nRemoved $removed duplicate sequences and wrote $file.nt.trim.unique.fasta\n";
	print $handle "Made mothur names file to keep track of duplicate sequences, $file.nt.trim.names\n\n";
} else {
	print $handle "\nUnique.seqs failed\nOperation aborted\n\n";
	print "\nUnique.seqs failed\nOperation aborted\n\n";
	die
}

# translates nucleotide data

my %fasta=readfasta("$file.nt.trim.unique.fasta",$handle);

my $fastaref=\%fasta;

translate($fastaref,"$file.aa.trim.unique.fasta",$handle);

if (-e "$file.aa.trim.unique.fasta") {
	print $handle "Translated nucleotide sequences and saved as $file.aa.trim.unique.fasta\n\n";
} else {
	print $handle "\nTranslation of nucleotide sequences failed\nOperation aborted\n\n";
	print "\nTranslation of nucleotide sequences failed\nOperation aborted\n\n";
	die
}

# runs amino acid sequences through muscle until file size stops changing OR starts looping between states

muscle("$file.aa.trim.unique.fasta","$file.aa.trim.unique.muscle.fasta",$muscleiterations,$handle,$logfile);

if (!-e "$file.aa.trim.unique.muscle.fasta") {
	print $handle "\nMUSCLE alignment failed\nOperation aborted\n\n";
	print "\nMUSCLE alignment failed\nOperation aborted\n\n";
	die
}

# back-translates aligned amino acids to nucleotides using original sequences

backtranslate("$file.aa.trim.unique.muscle.fasta","$file.nt.trim.unique.fasta","$file.nt.trim.unique.muscle.fasta");

if (-e "$file.nt.trim.unique.muscle.fasta") {
	print $handle "\nAmino acid alignment back-translated to nucleotides using $file.nt.trim.unique.fasta as a template.\n";
	print $handle "Writing back-translated alignment as $file.nt.trim.unique.muscle.fasta\n\n";
} else {
	print $handle "\nBack-translation failed\nOperation aborted\n\n";
	print "\nBack-translation failed\nOperation aborted\n\n";
	die
}

# If clustering is on, uses mothur to break retrieved and aligned sequences into OTUs, then prepares a taxonomy file with no ambiguities inside an OTU.

if ($cluster==1) {
	
	print $handle "\nBeginning mothur clustering operations\n\n";

	system("mothur \"#set.logfile(name=$logfile, append=T);dist.seqs(fasta=$file.nt.trim.unique.muscle.fasta, cutoff=".($cutoff*1.1).")\"");

	system("mothur \"#set.logfile(name=$logfile, append=T);cluster(column=$file.nt.trim.unique.muscle.dist, name=$file.nt.trim.names, hard=t, cutoff=$cutoff)\"");

	system("mothur \"#set.logfile(name=$logfile, append=T);get.oturep(fasta=$file.nt.trim.unique.muscle.fasta, column=$file.nt.trim.unique.muscle.dist, list=$file.nt.trim.unique.muscle.an.list, name=$file.nt.trim.names)\"");

	for(my $i = $cutoff;$i>=0.01;$i+=-0.01) {
		$i=sprintf("%0.2f",$i);
		unless (-e "$file.nt.trim.unique.muscle.an.$i.rep.fasta") {
			if ($i==0.01) {
				die "No appropriate OTU rep file\n";
			} else {next}
		}
		system("cp $file.nt.trim.unique.muscle.an.$i.rep.fasta $file.reference.fasta 1>>$logfile 2>>$logfile");
		system("cp $file.nt.trim.unique.muscle.an.$i.rep.names $file.reference.names 1>>$logfile 2>>$logfile");
		
		if (-e "$file.nt.trim.unique.muscle.dist") {print $handle "\nWrote distance matrix for $enzyme alignment as $file.nt.trim.unique.muscle.dist\n"}
		if (-e "$file.reference.fasta") {print $handle "\nWrote mothur reference alignment as $file.reference.fasta\n"}
		if (-e "$file.reference.names") {print $handle "\nWrote mothur reference names file as $file.reference.names\n"}
		else {
			print $handle "\nClustering operations failed\nOperation aborted\n\n";
			print "\nClustering operations failed\nOperation aborted\n\n";
			die
		}
			
 		system("rm $file.nt.trim.unique.muscle.an.* 1>>$logfile 2>>$logfile");	otutaxa($taxonomy,"$file.reference.names",'_',$file,"$file.reference.fasta","$file.template.fasta",$handle);
 		
		print $handle "\nEnding mothur clustering operations\n\n";
		last;
	}
} else {
# If clustering is off,
	print $handle "\nSkipping mothur clustering operations\n\n"; 
	system("cp $file.nt.trim.names $file.reference.names 1>>$logfile 2>>$logfile");
	system("cp $file.nt.trim.unique.muscle.fasta $file.reference.fasta 1>>$logfile 2>>$logfile");
	if (-e "$file.reference.fasta" && -e "$file.reference.names") {
		print $handle "\nWrote mothur reference alignment as $file.reference.fasta\n";
		print $handle "\nWrote mothur reference names file as $file.reference.names\n";
	} else {
		print $handle "\nFailed to write reference files\nOperation aborted\n\n";
		print "\nFailed to write reference files\nOperation aborted\n\n";
		die
	}
	otutaxa($taxonomy,"$file.reference.names",'_',$file,"$file.reference.fasta","$file.template.fasta",$handle);
	

}

system("makeblastdb -in $file.aa.trim.unique.fasta -dbtype prot -out $file.aa.db 1>>$logfile 2>>$logfile");

print $handle "\nWrote BLAST database $file.db\n\n";
print $handle "KEGGUR analysis complete\n";
close $handle;

exit;

#-----------------------------------------------------------------------------------------------------------
sub getlist 
{

# gets list of all annotated coding genes w. ARG0 handle from KEGG database and saves to ARG1
# Can use EC codes (give ec:code as ARG) or KEGG ortholog code (ko:code)
use strict;
use warnings;


my $code=$_[0];
my $file=$_[1];
my $logfile=$_[2];
my $count=0;

system('wget -O '.$file.'.temp.kegg http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+ko:'.$code." -a $logfile \n");

open(my $handle,"$file.temp.kegg");
	my @list=<$handle>;
close($handle);

system("rm $file.temp.kegg 1>>$logfile 2>>$logfile");

foreach my $element (@list) {
	$element =~ s/<[^>]*>//g;
        $element =~ s/([^ ][^ ]*)  *.*/$1/;
}
@list=grep(/^[a-z][a-z]*:[^ ]*/,@list);
@list=grep(/\S/,@list);

if ($file) {
open($handle,'>'.$file);
	foreach my $element (@list) {
		print $handle $element;
		$count ++;
	}
close($handle);
}

return $count;

}

#-------------------------------------------------------------------------------------------------------------
sub getkegg # gets nucleotide sequences of KEGG genes from list file ARG0 and saves them to ARG1
{

use strict;
use warnings;


my ($infile,$outfile,$logfile)=@_;
my @nt=();
my $count=0;

open(my $handle,$infile);
	my @list=<$handle>;
close($handle);

foreach my $element (@list) {
	chomp $element;
	system('wget -O '.$infile.'.temp.file http://www.genome.jp/dbget-bin/www_bget?-f+-n+n+'.$element." -a $logfile \n");
	open($handle,"$infile.temp.file");
	$count ++;
		my @ntseq=<$handle>;
	close($handle);
	system("rm $infile.temp.file 1>>$logfile 2>>$logfile");
	foreach my $seq (@ntseq) {
			$seq=~s/<[^>]*>//g;
			$seq=~s/^[^>atgcATGC].*//g;
			$seq=~s/\([^ ][^ ]*\)  *.*/$1/;
		}
	@ntseq=grep(/\S/,@ntseq);
	push(@nt,@ntseq);
}

if ($outfile) {
	open($handle,'>'.$outfile);
		foreach my $element (@nt) {
			print $handle $element;
		}
	close($handle);
}

return $count;

}
#----------------------------------------------------------------------------------------------------------
sub readfasta # return hash of headers and sequences from file ARG0
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
	    $sequence .= "\n";
	    if ($name) {
	    	$fasta{$name}=$sequence;
		$sequence='';
	    }
	    $name=$line;
	
       } else {
	    chomp($line);
            $sequence .= $line;
       }

}
$sequence .= "\n";

$fasta{$name}=$sequence;


return %fasta;
}

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
         print $handle "Bad codon \"$codon\"!!\n";
         return "X";
     }
 }

#-------------------------------------------------------------------------------------
sub translate {
# converts ARG0 fasta DNA hashref into peptide hash and writes to ARG1 file

use strict;
use warnings;


my ($ntfasta,$outfile,$handle)=@_;
my %aafasta=();

while (my ($name, $seq) = each %{$ntfasta}) {
	my $pep='';
	chomp($seq);
	my $inframe = length($seq);
	if (my $mod=$inframe % 3) {
		print $mod;
		$seq=substr($seq,0,$inframe-$mod);
	}
	for (my $ntp=0; $ntp < (length($seq)-2); $ntp=$ntp+3) {
                my $frame=substr($seq,$ntp,3);
                my $aa=codon2aa($frame,$handle);
            # If stop codon, don't include and stop translating
            if ($aa eq "*") {
            	last
            }
            else {
        		$pep=$pep.$aa;
        	}
        }
	$pep .= "\n";
	$aafasta{$name}=$pep;
}

if ($outfile) {
	my $href=\%aafasta;
	writefasta($href,$outfile);
}

return %aafasta;
}

#----------------------------------------------------------------------------------------
sub muscle {
# runs muscle on a set of sequences until the file sizes of 3 subsequent runs converge or ARG2 runs
# given fasta file as ARG0, outfile is ARG1
# requires MUSCLE to be in the system path

use strict;
use warnings;
use List::Util qw(min max first);


my ($infile,$outfile,$runs,$handle,$logfile)=@_;

my (@size,@loop,@loop2)=();
my $min=0;

system("muscle -in $infile -maxiters 9999 -out $infile.1.test.fasta -loga $logfile");
$size[0]=-s "$infile.1.test.fasta";
print $handle "\n\nMUSCLE iteration 1: file size $size[0]\n\n";

for(my $i=2; $i <= $runs; $i ++) {
	my $j=$i-1;
	system("muscle -in $infile.$j.test.fasta -maxiters 9999 -out $infile.$i.test.fasta -loga $logfile");
	$size[$j]=-s "$infile.$i.test.fasta";
	print $handle "\n\nMUSCLE iteration $i: file size $size[$j]\n\n";
	unless ($i<=2) {
		@loop=findloop(@size);
		if (@loop) {
			@loop2=splice(@loop,2,2+$loop[1]);
			$min=(min @loop2);
			my $idx = 1+$loop[0]+first{$loop2[$_] == $min} 0..$#loop2;
			my $targfile="$infile.$idx.test.fasta";
			system("mv $targfile $outfile 1>>$logfile 2>>$logfile");
			system("rm $infile.*.test.fasta 1>>$logfile 2>>$logfile");
			print $handle "\nAlignment converged in $idx iterations.\n";
			print $handle "MUSCLE iteration $idx chosen and written as $outfile\n\n";
			return;
		}
		else {print $handle "\nNo loop.\n\n"}
	}
}

print $handle "\nAlignment did not converge.\n";
print $handle "Final MUSCLE iteration ($runs) written as $outfile\n\n";
system("mv $infile.$runs.test.fasta $outfile 1>>$logfile 2>>$logfile");

return;

}

#-----------------------------------------------------------------------------------------
sub backtranslate {
# Take a translated amino acid fasta file and back-translates it to nucleotides, using a given
# nucleotide file. OK for AA file to be an alignment but not NT file.
# ARG[0]=amino acid filename, ARG[1]=nucleotide filename
# Returns results and also puts them into a file called ARG2.

use strict;
use warnings;


my ($aafile, $ntfile, $outfile)=@_;
my ($ntp, $aap)=0;
my (%aa, %nt, %bt)=();
my ($frame,$bttemp,$nttemp)='';

%aa=readfasta($aafile);
%nt=readfasta($ntfile);

while (my ($name,$seq) = each(%aa)) {
	chomp($seq);
	$nttemp=$nt{$name};
	for ($aap = 0; $aap <= (length($seq)-1); $aap++) {
		$frame=substr($seq,$aap,1);
		if ($frame =~ /-/) {
			$bttemp.="---";
			next;
		}
		elsif ($frame =~ /[A-Z*]/i) {
			$bttemp.=substr($nttemp,$ntp,3);
			$ntp+=3;
			next;
		}
		else {
			print "Bad Amino Acid\n";
			exit;
		}
	}
	$bt{$name}="$bttemp\n";
	$bttemp='';
	$ntp=0;

}

my $btref = \%bt;
writefasta($btref,$outfile);

return %bt;

}
#------------------------------------------------------------------
sub findloop {
# takes an array of values ARG0 and determines if there's a loop that propagates at least twice.
# Returns the parameters:
# mu = position where the loop starts (initial = pos 0)
# lambda = length of the loop
# the elements of the loop
# return variable @loop: 0 = mu, 1 = lambda, 2-2+lambda = elements

use strict;
use warnings;


my @loop = ();
my ($len,$confirm,$k,$flag,$lambda,$mu,$maxlam,$maxmu)=0;
my @size = @_;

$len=scalar @size;
$maxlam=int($len/3);

LOOPSEARCH:
for ($lambda=1; $lambda<=$maxlam; $lambda++) {
        $maxmu=$len+1-3*$lambda;
        for ($mu=0; $mu < $maxmu; $mu++) {
                if ($size[$mu]==$size[$mu+$lambda] && $size[$mu]==$size[$mu+2*$lambda]) {
                        $flag=1;
                        SUBLOOP:
                        for ($confirm=1; $confirm < $lambda; $confirm++) {
                                unless($size[$mu+$confirm]==$size[$mu+$lambda+$confirm] && $size[$mu+$confirm]==$size[$mu+2*$lambda+$confirm]) {
                                        $flag=0;
                                        last SUBLOOP
                                }
                        }
                if ($flag==1) {last LOOPSEARCH}
                } else {$flag=0}
        }
}

if ($flag==0) {
    return;
}

$loop[0]=$mu;
$loop[1]=$lambda;

for ($k = 2; $k < $lambda+2; $k++) {
        my $index=$mu+$k-2;
        $loop[$k]=$size[$index];
}

return @loop;

}

#----------------------------------------------

sub otutaxa {
# given a taxonomy file produced by getorgs.pl (arg 0) and a names file (arg 1),
# # checks to make sure taxonomies agree.  If they don't writes
# # a redundant taxonomy file that separates all found taxa by slashes
# # and a conservative taxonomy file that eliminates any ambiguous
# # taxa and only returns the highest-level unambiguous taxon
# # (gives "Unknown" if highest-level taxon is ambiguous)
# Arg2 is the separator character in the names file between the taxa code
# and the gene code/org name/whatever else
# Arg 3 is output file
# Arg 4 is a fasta file corresponding to the taxonomy file;
# If given, the program will remove all gaps to make a "template" file
# Arg 5 is fasta template output filename

use strict;
use warnings;
use List::Util qw(min max first);


my ($taxonomy,$namesfile,$sep,$outputfile,$fastain,$fastaout,$logfile)=@_;
my $maxtaxa=0;
my %listhash = ();
my (@taxonomy,@names,@count,@newtaxarray,@conservtaxarray)=();
my $handle;
my $redcount;

# writes gapless template fasta file

if (open($handle,$fastain)) {
	my @fasta = <$handle>;
	close($handle);
	open(my $handle,">$fastaout");
	foreach my $element (@fasta) {
		if ($element =~ /^>/) {print $handle $element}
		else {
			$element =~ s/[.-]//g;
			print $handle $element;
		}
	}
	print $logfile "\nWrote gapless mothur template file as $fastaout\n\n";
	close ($handle);
}

if (open ($handle,$taxonomy)) {
	@taxonomy = <$handle>;
	close ($handle);
} else {
	print "No such taxonomy file\n";
	print $logfile "Missing taxonomy file; run getorgs.pl\n\n";
	die
}

if (open ($handle,$namesfile)) {
	@names = <$handle>;
	close ($handle);
} else {
	print "No such names file\n";
	print $logfile "Missing names file.\n\n";
	die
}

foreach my $element (@taxonomy) {
        my $count = $element =~ tr/;//;
        push (@count,$count);
}

$maxtaxa=(max @count);
print $logfile "\nTaxonomy is $maxtaxa levels deep.\n";

# pads all taxonomy elements with na's out to maximum taxonomic depth
foreach my $element (@taxonomy) {
	chomp($element);
        my $count = $element =~ tr/;//;
        my $deficit = $maxtaxa-$count;
        $element .= ('na;' x $deficit);
}


# makes a hash of taxonomies keyed with organism names

my $aref=\@taxonomy;

my %tax = taxhash($aref);

# Gets all redundant sequences from names file and checks them using taxonomy hash to make sure their taxonomies agree
# Mushes up all taxonomic designations at each level for all redundant sequences in names file
# Then eliminates all that are identical
# Then makes two taxonomy files: one that keeps redundancies separated by slashes, another that deletes any taxonomic designations lower than or equal to the first redundant level.

foreach my $element (@names) {
	chomp($element);
        my $list = $element;
        $list =~ s/^.*\t//;
        my $name = $element;
	my $fullname = $name;
	$fullname =~ s/\t.*$//;
        $name =~ s/$sep.*$//;
        my @list = split(',',$list);
        foreach my $whoosit (@list) {
		$whoosit =~ s/$sep.*$//;
	}
        my $len = scalar @list;
        my %entryhash=();
        my @newtax = ();
        for (my $i=0 ; $i < $len; $i++) {
                my @lineage = split(';',$tax{$list[$i]});
                $entryhash{$i}=[@lineage];
        }

        for (my $i = 0; $i < $maxtaxa; $i++) {
                my @taxon=();
                for (my $j = 0; $j < $len; $j++) {
                        push (@taxon, ${$entryhash{$j}}[$i]);
                }
                @taxon = uniq(@taxon);
                my $taxon = join('/',@taxon);
                push (@newtax,"$taxon;");
        }
	my @trimtax=@newtax;
	while (pop(@trimtax) eq "na;") {pop(@newtax)}
        my $newtax = join('',@newtax);
        my $conservtax = $newtax;
        $conservtax =~ s/;[^;]*\/.*$/;/;
        $conservtax =~ s/;na;/;/g;
        if ($conservtax =~ /\//) {$conservtax = "Unknown;"}
        $newtax = "$fullname\t$newtax";
        $conservtax="$fullname\t$conservtax";
        if ($newtax ne $conservtax) {
        	$redcount += 1;
        	print $handle "\nDisagreement:\n$newtax\n$conservtax\n\n";
        }
	push (@newtaxarray,$newtax);
	push (@conservtaxarray,$conservtax);
}

if ($redcount) {
	print $logfile "\n$redcount taxonomies had at least 1 disagreement between redundant sequences.\n\n";
	open ($handle,">$outputfile.red.tax");
	foreach my $element (@newtaxarray) {print $handle "$element\n"}
	close ($handle);
	print $logfile "Redundant taxonomies saved as $outputfile.red.tax\n";
	
	open ($handle,">$outputfile.cons.tax");
	foreach my $element (@conservtaxarray) {print $handle "$element\n"}
	close ($handle);

	print $logfile "Conservative taxonomies with disagreements removed saved as $outputfile.cons.tax\n\n";
} else {
	print $logfile "\nNo taxonomy disagreements.\n\n";

	open ($handle,">$outputfile.cons.tax");
	foreach my $element (@conservtaxarray) {print $handle "$element\n"}
	close ($handle);

	print $logfile "Taxonomies saved as $outputfile.cons.tax\n\n";
}

return
}

sub uniq {
        return keys %{{ map { $_ => 1 } @_ }};
}

#------------------------------------------------------------------------------------

sub writefasta # takes referenced hash of headers and sequences ARG0 and writes a file or returns an array
{

use strict;
use warnings;


my ($fasta,$outfile,$logfile)=@_;
my @fasta = ();

while (my ($name,$seq) = each (%{$fasta})) {
	push(@fasta,$name);
	push(@fasta,$seq);
}

if ($outfile) {
	open (my $handle,'>'.$outfile);
		foreach my $element (@fasta) {
			print $handle $element;
		}
	close ($handle);
}

return @fasta;
}

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
	my $handle='';
	open ($handle,$input);
	@tax = <$handle>;
	close($handle);
} else {@tax = @{$input}}

foreach my $line (@tax) {
       my $code = $line;
       my $lineage = $line;
       $code =~ s/\|.*$//g;
       $lineage =~ s/^.*\t//g;
       $taxhash{$code}=$lineage;
}

return %taxhash;
}

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
