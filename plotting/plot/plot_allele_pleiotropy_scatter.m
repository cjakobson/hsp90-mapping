function [] = plot_allele_pleiotropy_scatter(dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


mapping_input=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);
mapping_input=mapping_input(mapping_input.isQtn==1,:);

unique_qtls=unique(mapping_input.index);

mapping_input.delta_deltaZ=mapping_input.deltaZbuffer-mapping_input.deltaZbaseline;

m=1;
for i=1:length(unique_qtls)
    
    temp_table=mapping_input(mapping_input.index==unique_qtls(i),:);
    
    for j=1:height(temp_table)
        
        for k=(j+1):height(temp_table)
            
            if (~strcmp(temp_table.condition{j},temp_table.condition{k}))
                v1(m)=temp_table.delta_deltaZ(j);
                v2(m)=temp_table.delta_deltaZ(k);
                m=m+1;
            end
            
        end
        
    end
    
end

hold on
scatter(v1,v2,10,'k','filled')
axis square
xlim([-0.5 0.5])
ylim([-0.5 0.5])
plot(xlim,[0 0],':r')
plot([0 0],ylim,':r')




end