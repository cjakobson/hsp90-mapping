function v_dist = calculate_telomere_distance_gene(dependency_directory,output_directory)


%background density of all genes
gene_data=readtable([dependency_directory 'TableS4.xls']);

telo_input=readtable([dependency_directory 'chromosome_length.txt']);


%convert chromosomes to numerals
chr_array={'I','II','III','IV','V','VI','VII','VIII','IX','X',...
    'XI','XII','XIII','XIV','XV','XVI'};


v_dist=nan(height(gene_data),1);

for i=1:length(v_dist)
   
    temp_chr_roman=gene_data.Chrom{i};
    temp_chr=find(ismember(chr_array,temp_chr_roman(4:end)));
    
    temp_pos=gene_data.SGD_Start(i);
    
    temp_chr_length=telo_input.Var3(telo_input.Var1==temp_chr);
    
    if ~isempty(temp_chr_length)
        temp_dist=abs(temp_chr_length-temp_pos);

        v_dist(i)=min([temp_pos,temp_dist]);
        %normalize to chromosome length
        v_dist(i)=v_dist(i)/temp_chr_length;
    end
    
end


end


