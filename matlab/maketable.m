function t  = maketable(in)
    
    gene_id = {};
    gene_symbol = {};
    num_not_SL = {};
    num_error = {};
    num_SL = {};
    
    
    for i = 1:length(in)
        gene_id{i} = num2str(in{i,1});
        gene_symbol{i} = in{i,2};
        num_not_SL{i} = in{i,3};
        num_error{i} = in{i,4};
        num_SL{i} = in{i,5};
    end
    
    gene_id = gene_id';
    gene_symbol = gene_symbol';
    num_not_SL = num_not_SL';
    num_error = num_error';
    num_SL = num_SL';
    
    t = table(gene_id, gene_symbol, num_not_SL, num_error, num_SL)
    
end
    
