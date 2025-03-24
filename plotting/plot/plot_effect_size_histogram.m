function [] = plot_effect_size_histogram(dependency_directory,output_directory)


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
to_plot{1}=abs(input_data_no_rad.deltaZbaseline);
to_plot{2}=abs(input_data_hsp90.deltaZbuffer-input_data_hsp90.deltaZbaseline);

%median(input_data_no_rad.varExp,'omitnan')
%median(input_data_hsp90.varExp,'omitnan')

histogram(to_plot{1},0:0.05:1,'Normalization','probability')
histogram(to_plot{2},0:0.05:1,'Normalization','probability')

ylim([0 0.3])
xlim([0 1])
legend({'no rad','Hsp90-dep.'})
ylabel('freq.')
xlabel('effect size')
axis square


end


