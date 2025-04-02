#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(shuffle);  # Import shuffle function

#For each species in your Go folder work out species name, then store in a hash with species name -> path to file
my @gos=`ls Go/*`;
my %go_key;
foreach my $sp (@gos){
    chomp $sp;
    my @split=split(/\./, $sp);
    my @sp_folder=split("\/", $split[0]);
    $go_key{$sp_folder[1]}=$sp;
}

# Ensure the correct number of arguments are provided
if (@ARGV != 3) {
    die "Usage: $0 <size_in_bp> <species><replicate numbers>\n";
}

my $size = $ARGV[0];
my $species = $ARGV[1];
my $replicates = $ARGV[2];

# Store the total scores and counts for each key
my %scores;

# Read all input files ending with '_gene_scores.txt'
my @files = glob("*.gene_scores.txt");

foreach my $file (@files) {
    open my $fh, '<', $file or die "Could not open file '$file' $!";
    print "Combining $file\n";
    while (<$fh>) {
        chomp;
        my ($col1, $col2, $score) = split /\t/;

        # Create a unique key for each row based on the first two columns
        my $key = "$col1\t$col2";

        # Add new score if less than prior:
        if ($scores{$key}){
            if ($scores{$key} eq "NA"){
                $scores{$key} = $score;
            }
            else{
                if ($scores{$key} >= $score){
                    $scores{$key} = $score;
                }
            }
            
        }
        else{
            $scores{$key} = $score;
        }
        
        
    }

    close $fh;
}

# Write the lowest scores to "Summary.tsv"
open my $out_fh, '>', "Summary.tsv" or die "Could not open file 'Summary.tsv' $!";

#Save output of genes tested:
my $back="$species\.$size\.bg.txt";
open my $backsave, '>', $back or die "Could not open file '$back': $!";

# Count of distance scores from break, e.g. <1 scores or <3 scores
my $distance_count=0;
my @closest_cutoff;
my @backgenes;
foreach my $key (keys %scores) {
    if ($scores{$key} eq "NA"){
        #Do nothing.
    }
    else{
        my @splh=split("\t", $key);
        if($splh[1] =~ m/\:/){
            my @sp1=split(/\:/, $splh[1]);
            $splh[1]=$sp1[1];
        }
        # Rename weird NCBI id rna- prefix
        if($splh[1] =~ m/rna-/){
                    my @sp1=split(/\-/, $splh[1]);
                    $splh[1]=$sp1[1];
                }
        if($splh[1] =~ m/-/){
            $splh[1] =~ s/\-/\_/g;
        }
        print $backsave "$splh[1]\n";
        push (@backgenes, "$splh[1]");
    }
}


for (my $i=0; $i<=$replicates; $i++){

    my $random_output_file = "$species.${size}\.random\.$i.txt";
    open my $closest_fh, '>', $random_output_file or die "Could not open file '$random_output_file': $!";

    # Shuffle the array and take the first $size elements
    my @random_genes = (shuffle @backgenes)[0 .. $size-1];

    print $closest_fh "$_\n" for @random_genes;

    close $closest_fh;

    print "Results random genes into $random_output_file\n";

    #Run GO enrichment analysis on lists

    print "ChopGO_VTS2_v12.pl -i $random_output_file --GO_file $go_key{$species} -bg $species\.$size\.bg.txt \n";
    `ChopGO_VTS2_v12.pl -i $random_output_file --GO_file $go_key{$species} -bg $species\.$size\.bg.txt`; 

}



