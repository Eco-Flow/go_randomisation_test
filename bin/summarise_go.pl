#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min);

# Parse command-line arguments
my $pvalue_type = "bonferroni";  # Default column
my $pvalue_cutoff = 0.05;        # Default threshold
GetOptions(
    "type=s"   => \$pvalue_type,
    "pvalue=f" => \$pvalue_cutoff
) or die "Usage: $0 --type <pvalue_column> --pvalue <cutoff>\n";

# Directory containing input files
my $dir = ".";  # Change if needed
my @files = glob("$dir/*.tab");  # Get all .tab files

# Extract the common prefix (species and number)
my $common_prefix;
if (@files) {
    ($common_prefix) = $files[0] =~ /^(.*?\.\d+)\./;  # Extract "Drosophila_santomea.200"
} else {
    die "No input files found.\n";
}

# Define output filenames using the extracted prefix
my $top50_output = "${common_prefix}.top50_lowest_${pvalue_type}_pvalues.txt";
my $summary_output = "${common_prefix}.go_term_occurrence_summary.txt";

open my $top50_fh, '>', $top50_output or die "Could not open '$top50_output': $!";
open my $summary_fh, '>', $summary_output or die "Could not open '$summary_output': $!";

# Hashes to store results
my %min_pvals;
my %term_count;
my %lowest_pvalues;
my %go_terms;  # Stores GO ID -> Term mapping

foreach my $file (@files) {
    open my $fh, '<', $file or die "Could not open '$file': $!";

    # Read header and find column index for the chosen p-value type
    my $header = <$fh>;
    chomp $header;
    my @columns = split(/\t/, $header);
    
    my %col_idx;
    for my $i (0 .. $#columns) {
        $col_idx{ lc $columns[$i] } = $i;  # Store index by lowercase column name
    }
    
    # Check if user-specified p-value type exists
    my $pval_col = $col_idx{ lc $pvalue_type };
    if (!defined $pval_col) {
        die "ERROR: Column '$pvalue_type' not found in header of $file\n";
    }

    while (my $line = <$fh>) {
        chomp $line;
        my @cols = split(/\t/, $line);

        # Extract relevant fields
        my $go_id = $cols[0];       # GO term ID
        my $term = $cols[1];        # Term description
        my $pvalue = $cols[$pval_col];  # Selected p-value column

        # Skip invalid lines
        next unless defined $pvalue and $pvalue =~ /^[\deE.-]+$/;

        # Convert scientific notation to decimal if necessary
        $pvalue = sprintf("%.10f", $pvalue) if $pvalue =~ /e/i;

        # Skip p-values above the user-defined cutoff
        next if $pvalue >= $pvalue_cutoff;

        # Store GO term name
        $go_terms{$go_id} = $term;

        # Update minimum p-value
        if (!exists $min_pvals{$go_id} or $pvalue < $min_pvals{$go_id}) {
            $min_pvals{$go_id} = $pvalue;
            $lowest_pvalues{$go_id} = [$term, $pvalue];  # Store term & p-value
        }

        # Count occurrences of this GO term
        $term_count{$go_id}++;
    }

    close $fh;
}

# Sort GO terms by minimum p-value and take top 50
my @top_50 = sort { $min_pvals{$a} <=> $min_pvals{$b} } keys %min_pvals;
@top_50 = @top_50[0..49] if @top_50 > 50;  # Keep only 50 terms

# Write Top 50 lowest p-values file
print $top50_fh "GO_ID\tTerm\t$pvalue_type\.p-value\n";
foreach my $go_id (@top_50) {
    my ($term, $pval) = @{$lowest_pvalues{$go_id}};
    print $top50_fh "$go_id\t$term\t$pval\n";
}

# Write GO term occurrence summary file
print $summary_fh "GO_ID\tTerm\tOccurrence_Count\n";
foreach my $go_id (sort { $term_count{$b} <=> $term_count{$a} } keys %term_count) {
    print $summary_fh "$go_id\t$go_terms{$go_id}\t$term_count{$go_id}\n";
}

close $top50_fh;
close $summary_fh;

print "Results saved to '$top50_output' and '$summary_output'\n";
