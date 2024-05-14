function [] = plot_cumulative_pval(dependency_directory,output_directory)



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

v_bins=0:50;

[~,v_type]=variant_types(input_data.variantType);

variant_labels={'protein-altering','synonymous','outside ORFs'};
mat_to_plot=nan(length(variant_labels),length(v_bins));

hold on
for i=1:length(variant_labels)
    
    variant_idx=v_type==i;
    
    temp_data=input_data(variant_idx,:);
    
    for j=1:length(v_bins)
        
        mat_to_plot(i,j)=sum(temp_data.pVal<=v_bins(j));
        
    end
    
    scatter(v_bins,mat_to_plot(i,:),25,'filled')
    
end

axis square
legend(variant_labels)
ylabel('cumulative freq.')
xlabel('p value')
xlim([0 max(v_bins)])
ylim([0 600])

end

