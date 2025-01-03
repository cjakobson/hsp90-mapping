function [] = plot_effect_size_abundance(dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

%just QTNs
input_data=input_data(input_data.isQtn==1,:);

[v_fraction_hsp90,v_type_hsp90]=variant_types(input_data.variantType);

v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);



abundance_data=readtable([dependency_directory '4932/4932-WHOLE_ORGANISM-integrated.txt']);

for i=1:height(abundance_data)
    
    temp_str=strsplit(abundance_data.x_string_external_id{i},'.');
    
    abundance_orf{i}=temp_str{2};
    
end

input_data.mean_abundance=nan(height(input_data),1);
for i=1:height(input_data)
    
    query_gene1=input_data.gene1{i};
    temp_abundance1=abundance_data.abundance(ismember(abundance_orf,query_gene1));
    
    if ~isempty(input_data.gene2{i})
    
        query_gene2=input_data.gene2{i};
        temp_abundance2=abundance_data.abundance(ismember(abundance_orf,query_gene2));
        
        input_data.mean_abundance(i)=mean([temp_abundance1,temp_abundance2]);
        
    else
        
        if ~isempty(temp_abundance1)
        
            input_data.mean_abundance(i)=temp_abundance1;
            
        end
        
    end
    
end

v1=input_data.mean_abundance;
v2=v_delta_delta_z;

scatter(v1,v2,10,'k','filled')
axis square
set(gca,'XScale','log')
ylim([0 1])
ylabel('|\Delta\DeltaZ|')
xlabel('est. protein abundance')
[r p]=corr(v1,v2,'rows','complete');
text(100,0.9,['r = ' num2str(r)])
text(100,0.8,['p = ' num2str(p)])


end

