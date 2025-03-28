function [] = plot_effect_size_ase_syn(dependency_directory,output_directory)


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

%only syn variants
[~,v_type]=variant_types(input_data.variantType);
input_data=input_data(v_type==2,:);

v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);

v_has_ase=logical(input_data.hasAseRad);

v_has_tag=logical(input_data.hasTagNoRad+input_data.hasTagRad);


to_plot{1}=v_delta_delta_z(logical(v_has_ase));
to_plot{2}=v_delta_delta_z(logical(v_has_tag.*~v_has_ase));


easy_box(to_plot)
ylim([0 0.5])
for i=1:length(to_plot)
    text(i,0.4,num2str(length(to_plot{i})))
end
for i=2:length(to_plot)
    [h p]=ttest2(to_plot{1},to_plot{i});
    text(i-0.5,0.45,num2str(p))
end
ylabel('\Delta\DeltaZ')
xticks(1:length(to_plot))
xtickangle(45)
xticklabels({'ASE','no ASE'})
title('synonymous variants')


end


