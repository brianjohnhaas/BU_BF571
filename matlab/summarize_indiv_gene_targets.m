% summarizeOtherTargets.m
% 
% Author: Iñigo Apaolaza
% Email: iaemparanza@ceit.es
% Date: 14/03/2017

close all
clear all
clc


top_files = textread('top_files.list', '%s');
n_top_files = length(top_files);



for ii = 1:n_top_files
    load(top_files{ii});
    clearvars -except ii TOP top_files n_GSM is_gMCS KO n_top_files
    n_KO = length(KO);
    n_is_gMCS = length(is_gMCS);
    fprintf('%d, GSM: %s, n_KO: %d, n_is_gMCS: %d\n', ii, top_files{ii}, n_KO, n_is_gMCS);
    TOP((ii-1)*n_KO+1:ii*n_KO, 1) = KO;
    TOP((ii-1)*n_KO+1:ii*n_KO, 2) = is_gMCS;
end


last_files = textread('last_files.list', '%s');
n_last_files = length(last_files);




for ii = 1:n_last_files
    load(last_files{ii})
    clearvars -except ii TOP LAST top_files last_files n_GSM is_gMCS KO n_top_files n_last_files
    n_KO = length(KO);
    n_is_gMCS = length(is_gMCS);
    fprintf('%d, GSM: %s, n_KO: %d, n_is_gMCS: %d\n', ii, last_files{ii}, n_KO, n_is_gMCS);
    LAST((ii-1)*n_KO+1:ii*n_KO, 1) = KO;
    LAST((ii-1)*n_KO+1:ii*n_KO, 2) = is_gMCS;
end

save('all_GSM_indiv_gene_results.mat', 'TOP', 'LAST', 'top_files', 'last_files')


clearvars -except TOP LAST

top_genes = unique(TOP(:, 1));
n_top_genes = length(top_genes);
last_genes = unique(LAST(:, 1));
n_last_genes = length(last_genes);

for i = 1:n_top_genes
    pos = find(TOP(:, 1) == top_genes(i));
    summary_TOP(i, 1) = top_genes(i);
    n_ok = sum(TOP(pos, 2) == 1);
    n_0 = sum(TOP(pos, 2) == 0);
    n_error = sum(TOP(pos, 2) == -1);
    summary_TOP(i, 2) = n_ok;
    summary_TOP(i, 3) = n_error;
    summary_TOP(i, 4) = n_0;
end
pos = summary_TOP(:, 2) == 0 & summary_TOP(:, 3) == 0;
summary_TOP = [summary_TOP(~pos, :); summary_TOP(pos, :)];

for i = 1:n_last_genes
    pos = find(LAST(:, 1) == last_genes(i));
    summary_LAST(i, 1) = last_genes(i);
    n_ok = sum(LAST(pos, 2) == 1);
    n_0 = sum(LAST(pos, 2) == 0);
    n_error = sum(LAST(pos, 2) == -1);
    summary_LAST(i, 2) = n_ok;
    summary_LAST(i, 3) = n_error;
    summary_LAST(i, 4) = n_0;
end
pos = summary_LAST(:, 2) == 0 & summary_LAST(:, 3) == 0;
summary_LAST = [summary_LAST(~pos, :); summary_LAST(pos, :)];
clear TOP LAST

load('symbol_vs_entrez.mat');
[~, ~, ind] = intersect(summary_TOP(:, 1), symbol_entrez.entrez, 'stable');
TOP = [mat2cell(summary_TOP(:, 1), ones(n_top_genes, 1), 1), symbol_entrez.symbol(ind), mat2cell(summary_TOP(:, 2:end), ones(n_top_genes, 1), ones(3, 1))];
[~, ~, ind] = intersect(summary_LAST(:, 1), symbol_entrez.entrez, 'stable');
LAST = [mat2cell(summary_LAST(:, 1), ones(n_last_genes, 1), 1), symbol_entrez.symbol(ind), mat2cell(summary_LAST(:, 2:end), ones(n_last_genes, 1), ones(3, 1))];

save('indiv_results.mat', 'TOP', 'LAST')

% xlswrite('summary_otherTargets_30clines.xlsx', TOP, 'Achilles Essential');
% xlswrite('summary_otherTargets_30clines.xlsx', LAST, 'Achilles Not Essential');
