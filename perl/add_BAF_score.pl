#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


my $usage = "\n\n\tusage: $0 top_or_last_genes.csv.wSym\n\n";

my $infile = $ARGV[0] or die $usage;

main: {

    print STDERR "-getting sample name mappings\n";
    my %GSM_to_sample_name = &parse_sample_name_mappings();

    print STDERR "-parsing BAF score matrix\n";
    my %BAF_scores = &parse_BAF_scores();

    print STDERR "-processing outputs\n";
    open(my $fh, $infile) or die "Error, cannot open file: $infile";
    my $header = <$fh>;
    chomp $header;
    print "sample_name\t$header\tBAF\n";
    while(<$fh>) {
        chomp;
        my ($gsm, $gene_sym, $gene_id, $cplex) = split(/\t/);
        my $sample_name = $GSM_to_sample_name{$gsm} or die "Error, no sample name for $gsm";
        my $baf = $BAF_scores{$sample_name}->{$gene_sym};
        if (! defined $baf) {
            $baf = "NA";
        }

        print join("\t", $sample_name, $gsm, $gene_sym ,$gene_id, $cplex, $baf) . "\n";

    }
    close $fh;


    exit(0);
    
}

####
sub parse_BAF_scores {
    my $matrix = "$FindBin::Bin/../gene_BayesFactor_scores.matrix";

    open(my $fh, $matrix) or die $!;
    my $header = <$fh>;
    chomp $header;
    my @cell_lines = split(/\t/, $header);
        
    my %sample_gene_to_BAF;

    while(<$fh>) {
        chomp;
        my @x = split("\t");
        my $gene_sym = shift @x;

        for (my $i = 0; $i <= $#x; $i++) {
            my $cell_line = $cell_lines[$i];
            my $val = $x[$i];
            $sample_gene_to_BAF{$cell_line}->{$gene_sym} = $val;
        }
    }
    close $fh;

    return(%sample_gene_to_BAF);
}
    

####
sub parse_sample_name_mappings {
    

    my %mappings;

    my $sample_name_mappings_file = "$FindBin::Bin/../CCLE/sample_name_mappings.txt";
    open(my $fh, $sample_name_mappings_file) or die $!;
    while(<$fh>) {
        chomp;
        my ($cell_line, $gsm) = split("\t");
        $mappings{$gsm} = $cell_line;
    }
    close $fh;

    return(%mappings);
}
