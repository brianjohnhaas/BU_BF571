% summarizeOtherTargets.m
% 
% Author: Iñigo Apaolaza
% Email: iaemparanza@ceit.es
% Date: 14/03/2017

close all
clear all
clc

GSM_list{1, 1} = 'GSM886864';
GSM_list{2, 1} = 'GSM886891';
GSM_list{3, 1} = 'GSM886892';
GSM_list{4, 1} = 'GSM886963';
GSM_list{5, 1} = 'GSM886997';
GSM_list{6, 1} = 'GSM886999';
GSM_list{7, 1} = 'GSM887000';
GSM_list{8, 1} = 'GSM887027';
GSM_list{9, 1} = 'GSM887046';
GSM_list{10, 1} = 'GSM887049';
GSM_list{11, 1} = 'GSM887058';
GSM_list{12, 1} = 'GSM887141';
GSM_list{13, 1} = 'GSM887142';
GSM_list{14, 1} = 'GSM887175';
GSM_list{15, 1} = 'GSM887262';
GSM_list{16, 1} = 'GSM887291';
GSM_list{17, 1} = 'GSM887300';
GSM_list{18, 1} = 'GSM887320';
GSM_list{19, 1} = 'GSM887338';
GSM_list{20, 1} = 'GSM887355';
GSM_list{21, 1} = 'GSM887364';
GSM_list{22, 1} = 'GSM887421';
GSM_list{23, 1} = 'GSM887441';
GSM_list{24, 1} = 'GSM887496';
GSM_list{25, 1} = 'GSM887499';
GSM_list{26, 1} = 'GSM887541';
GSM_list{27, 1} = 'GSM887576';
GSM_list{28, 1} = 'GSM887640';
GSM_list{29, 1} = 'GSM887691';
GSM_list{30, 1} = 'GSM887751';
n_GSM = length(GSM_list);

for ii = 1:n_GSM
    load(fullfile('.', 'data', 'output', ['BF_isgMCS_' GSM_list{ii} '_TOP_5min.mat']));
    clearvars -except ii TOP GSM_list n_GSM is_gMCS KO TOP_GSM
    n_KO = length(KO);
    n_is_gMCS = length(is_gMCS);
    fprintf('%d, GSM: %s, n_KO: %d, n_is_gMCS: %d\n', ii, GSM_list{ii}, n_KO, n_is_gMCS);
    TOP((ii-1)*n_KO+1:ii*n_KO, 1) = KO;
    TOP((ii-1)*n_KO+1:ii*n_KO, 2) = is_gMCS;

    for j = (ii-1)*n_KO+1:ii*n_KO
                TOP_GSM(j) = {GSM_list{ii,1}};
    end

end

for ii = 1:n_GSM
    load(fullfile('.', 'data', 'output', ['BF_isgMCS_' GSM_list{ii} '_LAST_5min.mat']));
    clearvars -except ii TOP LAST GSM_list n_GSM is_gMCS KO TOP_GSM LAST_GSM
    n_KO = length(KO);
    n_is_gMCS = length(is_gMCS);
    fprintf('%d, GSM: %s, n_KO: %d, n_is_gMCS: %d\n', ii, GSM_list{ii}, n_KO, n_is_gMCS);
    LAST((ii-1)*n_KO+1:ii*n_KO, 1) = KO;
    LAST((ii-1)*n_KO+1:ii*n_KO, 2) = is_gMCS;

    for j = (ii-1)*n_KO+1:ii*n_KO
                LAST_GSM(j) = {GSM_list{ii,1}};
    end


end

save('all_GSM_gene_results.mat', 'TOP', 'LAST', 'TOP_GSM', 'LAST_GSM');

    
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

save('results.mat', 'TOP', 'LAST')


% xlswrite('summary_otherTargets_30clines.xlsx', TOP, 'Achilles Essential');
% xlswrite('summary_otherTargets_30clines.xlsx', LAST, 'Achilles Not Essential');
