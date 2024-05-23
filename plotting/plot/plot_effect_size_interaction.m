function [] = plot_effect_size_interaction(dependency_directory,output_directory)


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

v_buffer=abs(input_data.deltaZbuffer)>=abs(input_data.deltaZbaseline);

v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);

to_plot{1}=v_delta_delta_z(v_buffer);
to_plot{2}=v_delta_delta_z(~v_buffer);

hold on
easy_box(to_plot)
xlim([0.5 2.5])
ylim([0 0.5])
[h p]=ttest2(to_plot{1},to_plot{2});
text(1.5,0.4,num2str(p))
ylabel('\Delta\DeltaZ')
for i=1:length(to_plot)
    text(i,0.45,num2str(length(to_plot{i})))
end


end


