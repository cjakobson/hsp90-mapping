function [] = plot_density_telomere_type(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

v_dist = calculate_telomere_distance(dependency_directory,output_directory);

%subset by coding vs regulatory
[~,v_type]=variant_types(variant_info.variantType);

%just QTNs
input_data=input_data(input_data.isQtn==1,:);
[~,v_type_mapping]=variant_types(input_data.variantType);


v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);

v_buffer=abs(input_data.deltaZbuffer)>=abs(input_data.deltaZbaseline);

v_dist_mapping=v_dist(input_data.index);

%just plot regulatory
v_dist=v_dist(v_type==3);
v_dist_mapping=v_dist_mapping(v_type_mapping==3);


%also plot density
v1=v_dist;
v2=v_dist_mapping;

v_hist_bins=0:0.05:0.5;

hold on
histogram(v1,v_hist_bins,'Normalization','probability')
histogram(v2,v_hist_bins,'Normalization','probability')
axis square
xlim([0 0.5])
ylim([0 0.25])
ylabel('frequency')
[h p]=kstest2(v1,v2);
text(0.4,0.1,num2str(p))


end


