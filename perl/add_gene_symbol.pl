#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;


my $usage = "usage: $0 genes.csv";

my $genes_csv = $ARGV[0] or die $usage;

main: {

    my %gene_id_to_symbol = &get_gene_symbol_info();

    
    open(my $fh, $genes_csv) or die "Error, cannot open file: $genes_csv";
    my $header = <$fh>;
    print join("\t", "sample", "sym", "id", "cplex") . "\n";
    
    while (<$fh>) {
        chomp;
        my ($gsm, $id, $cplex) = split(/,/);
        my $gene_sym = $gene_id_to_symbol{$id} or die "Error, no gene symbol for $id";

        print join("\t", $gsm, $gene_sym, $id, $cplex) . "\n";
        
    }

    exit(0);
}

####
sub get_gene_symbol_info {

    my %info;
    
    my $gene_symbol_info_file = "$FindBin::Bin/../symbol_vs_entrez.tsv";
    open(my $fh, $gene_symbol_info_file) or die "Error, cannot open file: $gene_symbol_info_file";
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($id, $symbol) = @x;
        $info{$id} = $symbol;
    }
    close $fh;

    return(%info);
}
