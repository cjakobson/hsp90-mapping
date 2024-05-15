function [] = plot_effect_size_line_crossing(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



input_data_hsp90=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);
input_data_hsp90=input_data_hsp90(input_data_hsp90.isQtn==1,:);

input_data_no_rad=readtable([dependency_directory 'linear_no_rad_fdr_0.05.csv']);
input_data_no_rad=input_data_no_rad(input_data_no_rad.isQtn==1,:);

hold on
histogram(abs(input_data_no_rad.deltaZbaseline),0:0.05:1,'Normalization','probability')
histogram(abs(input_data_hsp90.deltaZbuffer-input_data_hsp90.deltaZbaseline),0:0.05:1,'Normalization','probability')

ylim([0 0.3])
xlim([0 1])
legend({'no rad','Hsp90-dep.'})


end


