function output_table = parse_tss_tes(input_table,dependency_directory,output_directory)

tss_data=readtable([dependency_directory 'TableS4.xls']);

[~,v_type]=variant_types(input_table.variantType);
v_type=v_type';

input_table.tssDistGene1=nan(height(input_table),1);
input_table.tssDistGene2=nan(height(input_table),1);
input_table.tesDistGene1=nan(height(input_table),1);
input_table.tesDistGene2=nan(height(input_table),1);

input_table.upstreamGene1=nan(height(input_table),1);
input_table.upstreamGene2=nan(height(input_table),1);

for i=1:height(input_table)
    
    if v_type(i)==3
        
        query_gene1=input_table.gene1{i};
        query_gene2=input_table.gene2{i};
        
        gene_idx1=find(ismember(tss_data.Name,query_gene1));
        gene_idx2=find(ismember(tss_data.Name,query_gene2));
        
        if ~isempty(gene_idx1)
            
            temp_tss1=tss_data.x5__UTR_Start(gene_idx1);
            temp_tss_dist1=input_table.pos(i)-temp_tss1;
            
            temp_tes1=tss_data.x3__UTR_End(gene_idx1);
            temp_tes_dist1=input_table.pos(i)-temp_tes1;
            
            if length(query_gene1)>=7
                
                if strcmp(query_gene1(7),'C')
                    
                    temp_tss_dist1=-temp_tss_dist1;
                    temp_tes_dist1=-temp_tes_dist1;
                    
                end
                
            end
            
            input_table.tssDistGene1(i)=temp_tss_dist1;
            input_table.tesDistGene1(i)=temp_tes_dist1;
            
        end
        
        if ~isempty(gene_idx2)
            
            temp_tss2=tss_data.x5__UTR_Start(gene_idx2);
            temp_tss_dist2=input_table.pos(i)-temp_tss2;
            
            temp_tes2=tss_data.x3__UTR_End(gene_idx2);
            temp_tes_dist2=input_table.pos(i)-temp_tes2;
            
            if length(query_gene2)>=7
                
                if strcmp(query_gene2(7),'C')
                    
                    temp_tss_dist2=-temp_tss_dist2;
                    temp_tes_dist2=-temp_tes_dist2;
                    
                end
                
            end
            
            input_table.tssDistGene2(i)=temp_tss_dist2;
            input_table.tesDistGene2(i)=temp_tes_dist2;
            
        end
        
        if abs(input_table.tssDistGene1(i))<abs(input_table.tesDistGene1(i))
            
            input_table.upstreamGene1(i)=1;
            
        else
            
            input_table.upstreamGene1(i)=0;
            
        end
        
        if abs(input_table.tssDistGene2(i))<abs(input_table.tesDistGene2(i))
            
            input_table.upstreamGene2(i)=1;
            
        else
            
            input_table.upstreamGene2(i)=0;
            
        end
        
    end
    
end

output_table=input_table;


end