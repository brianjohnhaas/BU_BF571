#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use FindBin;
use Carp;
use lib ($ENV{EUK_MODULES});
use DelimParser;

my $usage = "usage: $0\n\n";

my $data_home = "$FindBin::Bin";


my $matlab_script = &get_matlab_script();

sub main {
    
    system("mkdir -p data/output");
    system("mkdir -p data/input/G_matrices");

    open(my $fh, "new_top_last_candidates_via_BAFnExpr.tsv") or die $!;
    my $delim_parser = new DelimParser::Reader($fh, "\t");

    my %data;
    
    my %gsm_to_top_or_last_KOs;
    while (my $row = $delim_parser->get_row() ) {
        my $top_or_last = $row->{'newTOP|newLAST'};
        my $gene_id = $row->{gene_id};
        my $gsm = $row->{GSM};

        push (@{$data{$gsm}->{$top_or_last}}, $gene_id);
    }


    foreach my $gsm (keys %data) {
        foreach my $top_or_last (keys %{$data{$gsm}}) {

            my @KOs = @{$data{$gsm}->{$top_or_last}};
            my $KO_list = join("; ", @KOs);
            my $GSM_name = $gsm;
            
            my $script = $matlab_script;
            $script =~ s/__KO_LISTING__/$KO_list/;
            $script =~ s/__GSM_CSV_FILENAME__/$GSM_name.csv/;
            $script =~ s/__TOP_OR_LAST__/$top_or_last/;
            
            my $scriptname = "${GSM_name}_$top_or_last.m";
            open(my $ofh, ">$scriptname") or confess "Error, cannot write to file: $scriptname";
            print $ofh $script;
            close $ofh;
            print STDERR "-wrote: $scriptname\n";
        }
        
    }
    close $fh;
        
    # use  .matlab-2012a
    # matlab -nodesktop -nodisplay -nosplash -r "run('/seq/RNASEQ/bf5seg/small_test/test.m');quit"
    
        
    exit(0);
}



sub get_matlab_script {

my $matlab_script = qq^

close all;
clear all;
clc;

% Add paths
addpath('/broad/software/free/Linux/redhat_6_x86_64/pkgs/cplex_12.6/cplex/matlab/');
addpath('$data_home/matlablib/');
warning('off');

% Inputs
GSM_name =  '__GSM_CSV_FILENAME__';
KO = [__KO_LISTING__];
% example:  280; 31; 205; 64499; 54575; 54600; 2678; 7177; 2746; 276; 277; 278; 6817; 1608; 279; 54659; 54657; 56052; 54576; 54579
n_KO = length(KO);
target_b = 1e-3;
timelimit = 5*60;

% Load Network
load('$data_home/Recon2v04_RPMI1640_purgedforMCS.mat');
nbio = find(strcmp(Network.rxndata.name, 'Generic human biomass reaction'));
[n_mets, n_rxns] = size(Network.S);

% Load GSMs
GSM = importdata(fullfile('$data_home/fRMA_summarized', GSM_name));
expression_binary_1(:, 1) = GSM.data(:, 1);
z = GSM.data(:, 2);
expression_binary_1(:, 2) = double(z > 5);
expression_binary_1(expression_binary_1(:, 2) == 0, 2) = -1;

% Get the ENTREZ ID of the genes in Recon2
ENTREZ_ID = Network.genes;
n_ENTREZ_ID = length(ENTREZ_ID);
for i = 1:n_ENTREZ_ID
    pos_dot = ismember(ENTREZ_ID{i}, '.');
    pos_dot = find(pos_dot);
    ENTREZ_ID{i} = ENTREZ_ID{i}(1:pos_dot-1);
end
ENTREZ_ID = cellfun(\@str2num, ENTREZ_ID);
ENTREZ_ID = unique(ENTREZ_ID);
ENTREZ_ID = sort(ENTREZ_ID, 'ascend');

% Keep only the expression_binary of the genes in Recon2
[genes, i_a, i_b] = intersect(ENTREZ_ID, expression_binary_1(:, 1));
expression_binary(:, 1) = ENTREZ_ID;
expression_binary(i_a, 2) = expression_binary_1(i_b, 2);

% Change the genes in the M set to the H set
expression_binary(expression_binary(:, 2) == 0, 2) = 1;

for i = 1:n_KO
    act_expression_binary = expression_binary;
    act_KO = KO(i);
    G_file = ['G_' GSM_name(1:end-4) '_KO-' num2str(act_KO) '.mat'];
    % Set the expression of the KO gene to -1 in order to enable its knock-out
    act_expression_binary(act_expression_binary(:, 1) == act_KO, 2) = -1;

    % Compute MCSs
    [is_gMCS(i, 1), cplex{i, 1}] = isGeneMCS(Network, nbio, target_b, expression_binary, act_KO, timelimit, G_file);
    save('-v7.3', fullfile('.', 'data', 'output', ['BF_isgMCS_' GSM_name(1:end-4) '___TOP_OR_LAST___' num2str(timelimit/60) 'min.mat']));
end

^;

    return($matlab_script);
}


&main();
