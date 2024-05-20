function [] = plot_gene_age_interaction(buffer_switch,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

input_data=input_data(input_data.isQtn==1,:);

v_buffer=abs(input_data.deltaZbuffer)>=abs(input_data.deltaZbaseline);

if buffer_switch==1
    input_data=input_data(v_buffer,:);
elseif buffer_switch==2
    input_data=input_data(~v_buffer,:);
end


v_delta_deltaZ=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);

group_labels={'Group I','Group II','Group III','WGD','Group IV','Group V'};

%classify interactions
for i=1:length(group_labels)
    
    v_temp=logical((input_data.gene1age==i)+(input_data.gene2age==i));
    
    to_plot{i}=v_delta_deltaZ(v_temp);
    
end

hold on
easy_box(to_plot)
set(gca,'YScale','log')
ylim([1e-2 1])
xticks(1:length(group_labels))
xtickangle(45)
xticklabels(group_labels)
for i=1:length(group_labels)
    
    text(i,0.8,num2str(length(to_plot{i})))
    [p h]=ranksum(to_plot{i},to_plot{end});
    text(i,0.6,num2str(p))
    
end
axis square
  


end


