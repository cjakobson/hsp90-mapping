function [] = plot_gene_density_telomere(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

v_dist = calculate_telomere_distance(dependency_directory,output_directory);

%just QTNs
input_data=input_data(input_data.isQtn==1,:);

v_dist_mapping=v_dist(input_data.index);


v_dist_gene= calculate_telomere_distance_gene(dependency_directory,output_directory);


%also plot density
v1=v_dist_gene;
v2=v_dist_mapping;

v_hist_bins=0:0.05:0.5;

hold on
histogram(v1,v_hist_bins,'Normalization','probability')
histogram(v2,v_hist_bins,'Normalization','probability')
axis square
xlim([0 0.5])
ylim([0 0.2])
ylabel('frequency')
[h p]=kstest2(v1,v2);
text(0.4,0.1,num2str(p))


end


