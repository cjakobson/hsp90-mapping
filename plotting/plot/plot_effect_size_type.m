function [] = plot_effect_size_type(dependency_directory,output_directory)


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

v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);

[~,v_type]=variant_types(input_data.variantType);

variant_labels={'protein-altering','synonymous','outside ORFs'};

for i=1:length(variant_labels)

    to_plot{i}=v_delta_delta_z(v_type==i);

end


hold on
easy_box(to_plot)
xlim([0.5 3.5])
ylim([0 0.5])
for i=1:length(to_plot)
    for j=(i+1):length(to_plot)
        [h p]=ttest2(to_plot{i},to_plot{j});
        text((i+j)/2,0.3+0.1*i,num2str(p))
    end
end
ylabel('\Delta\DeltaZ')


end


