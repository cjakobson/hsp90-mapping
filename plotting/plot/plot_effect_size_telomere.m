function [] = plot_effect_size_telomere(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

v_dist = calculate_telomere_distance(dependency_directory,output_directory);

v_dist=v_dist;

%just QTNs
input_data=input_data(input_data.isQtn==1,:);

v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);

v_buffer=abs(input_data.deltaZbuffer)>=abs(input_data.deltaZbaseline);

v_dist_mapping=v_dist(input_data.index);


%binned boxplot by interaction type
%v_bins=0:1e5:5e5;
%v_bins=[0 5e4 2e5 1e6];
v_bins=[0 0.05 0.25 0.5];

v1=v_dist_mapping;
v2=v_delta_delta_z;

m=1;
for i=1:(length(v_bins)-1)
    
    temp_idx=logical((v1>v_bins(i)).*(v1<=v_bins(i+1)));
    
    to_plot{m}=v2(logical(temp_idx.*v_buffer));
    m=m+1;
    
    to_plot{m}=v2(logical(temp_idx.*~v_buffer));
    m=m+1;
    
end


easy_box(to_plot)
ylim([1e-2 1])
set(gca,'YScale','log')
axis square

for i=1:(length(v_bins)-1)
    [p h]=ranksum(to_plot{2*(i-1)+1},to_plot{2*(i-1)+2});
    text(2*(i-1)+1.2,0.8,num2str(p))
end
temp_labels=cellfun(@length,to_plot);
for i=1:length(temp_labels)
    text(i,0.02,num2str(temp_labels(i)))
end
ylabel('|\Delta\DeltaZ|')



end


