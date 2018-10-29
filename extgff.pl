#!/usr/bin/perl
####################################################################################################
### Extracting GFF Features from Sequences                                                       ###
### $sh .pl drome_genes.gff drome_seq.fa                                                         ###
###                                                                                              ###
### This program retrieve those substrings from a DNA sequence that are delimited                ###
### by the start and end coordinates of a GFF record.                                            ###
###                                                                                              ###
### It also calculates the GC content of each sequence and total GC content                      ###
###                                                                                              ###
### It returns a tab delimited output:                                                           ###
### ID.group	start_position	end_position	GC percentage	Sequence			 ###
###                                                                                              ###
### 				Pozo Ocampo, Fernando                       			 ###                	
### 										December 5, 2016 ###           	
####################################################################################################

use strict;
use warnings;

&arguments_number;
&read_gff;
my $finalsequence = &read_fasta; # Getting last value from subroutine - > $dna = whole sequence without header and blank spaces in one line

my ($n, $f, $e, $dnaoutput, $gcpercentage, $totalgcpercentage, $dnareverse); 
my @names, my @first, my @final, my @strands;
my $identifiers = scalar (@names);  # @names has to be interpreted in scalar context in my loop

for ($n = 0; $n < $identifiers; $n++) {
	
    	$f = ($first [$n]) - 1; 		# Sequence correction 1
    	$e = ($final[$n] - $first[$n]) + 1; 	# Sequence correction 2
	$dnaoutput = substr ($finalsequence, $f, $e);	# Substring -> parsing the genome and taking nucleotides from each group position. It will be performed for + and -
	$gcpercentage = &gc_counter($dnaoutput);       # Counting the percentage of GC nucleotides. It will be performed for + and -. GC content is the same if you make revcom or not

	# if ($names[$n] =~ /^exon/) {  # If you want to print only one or more ID feature...

    	if ($strands[$n] eq "+") {	# If our array @strands is the direct chain...
		print STDOUT "$names[$n]\t$first[$n]\t$final[$n]\tGC-content: $gcpercentage%\t$dnaoutput\n";
    	} # printing the results of my arrays and my 2 strings
    	elsif ($strands[$n] eq "-"){
	$dnareverse = &revcom ($dnaoutput); # Reversing and getting the complementary chain of your DNA
		print STDOUT "$names[$n]\t$first[$n]\t$final[$n]\tGC-content: $gcpercentage%\t$dnareverse\n";
    	}
	# } Features curly bracket
};

$totalgcpercentage = &gc_counter($finalsequence);	# Calling the subroutine and printing GC of all genome
print "\n******The whole sequence have $totalgcpercentage % of GC content******\n";


sub arguments_number {		# Check if the number of arguments is 2 : file gff, file fasta
	if (@ARGV > 2)  {
	die "Too much arguments\n";
	}
	elsif (@ARGV == 1) { 
	die "Not enought arguments\n";
	}
};


sub read_fasta {

	my $fasta_file = $ARGV[1];
	my ($sequence, $dna);

open(FASTA, $fasta_file) || die "Can't open fasta file $!";
	while (<FASTA>) {

  		next if /^>/o;	 
    		chomp;
    		$sequence = substr($_, 0, length($_));

    		unless (defined $dna){	# Introducing the sequence in string $dna...
		$dna = $sequence;		 
    		}
    		else {
		$dna .= $sequence;
    		}
	}
close (FASTA);

return $dna;
};


sub read_gff {			# Initilizing a subroutine, it opens a gff file and take into an array the columns of interest

	my $gff_file = $ARGV[0]; 
	my ($name, $type, $start, $end, $strand);

open(GFF, $gff_file) || die "Can't open gff file $!";
	while (<GFF>) {

   	next if /^\s*$/o;	# Starting the next iteration of the loop 
    	chomp; 			# Checking empty lines

    	$name =  (split(/\s/, $_))[2];	# Spliting each column aiming to take the position of created array
    	$type = (split(/\s/, $_))[9];
    	$type =~ tr/";//d;		# Removing " and ; from ID
	$name .= ".$type";
    	push(@names, $name); # Adding two values (exon.group) to the first array

    	$start = (split(/\s/, $_))[3];
    	push(@first, $start);		# Adding start position to the second array 

    	$end = (split(/\s/, $_))[4];
    	push(@final, $end);

    	$strand = (split(/\s/, $_))[6];	# Adding strand position to the last array created
    	push(@strands, $strand);
	}
close (GFF);		# These variables have not been implement in a local enviroment
};

sub revcom {		 
	my($dnago) = $_[0];	# It contains the parameter passed to this subroutine
	my $revcom  = reverse $dnago; # First reverse the sequence
																					$revcom =~ tr/ACGTacgt/TGCAtgca/; # Next, complement the sequence, dealing with upper and lower case # A->T, T->A, C->G, G->C
return $revcom;			# Returning variable from this subroutine
};


sub gc_counter {

	my ($dnacount) = $_[0];	
	my ($char, $C, $G, $GC, $gc100);
	my $seqlen = length($dnacount);

$C = $G = 0;			# Important !! It has to start the counting by 0

for (my $n = -1; $n < $seqlen; $n++) {	#### WARNING: Sequence starts in -1 position of the string. It changes from ($n = 0 ?). It 						performs in a right way in order to count properly each sequence 
 	$char = uc(substr($dnacount, $n, 1));
   	SWITCH: {
		$char eq 'G' && ($G++, last SWITCH);
		$char eq 'C' && ($C++, last SWITCH);
	$GC = ($G + $C); 			# Total sum of my count
	$gc100 = ($GC / ($seqlen)) * 100;	# Calculating the percentage 
		}
	}
return $gc100;		
# return $GC;		# Optional returning variable from this subroutine (GC nucleotides). Better if you choose only one of them
};
