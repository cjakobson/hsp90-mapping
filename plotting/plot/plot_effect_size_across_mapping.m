function [] = plot_effect_size_across_mapping(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

hsp90_input=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);
no_rad_input=readtable([dependency_directory 'linear_no_rad_fdr_0.05.csv']);
rad_input=readtable([dependency_directory 'linear_rad_fdr_0.05.csv']);

hsp90_qtn_idx=hsp90_input.isQtn==1;
v_buffer=abs(hsp90_input.deltaZbuffer)>=abs(hsp90_input.deltaZbaseline);
v_potentiate=abs(hsp90_input.deltaZbuffer)<abs(hsp90_input.deltaZbaseline);

clear to_plot
to_plot{1}=abs(hsp90_input.deltaZbaseline(logical(hsp90_qtn_idx.*v_buffer)));
to_plot{2}=abs(hsp90_input.deltaZbaseline(logical(hsp90_qtn_idx.*v_potentiate)));

to_plot{3}=abs(no_rad_input.deltaZbaseline(no_rad_input.isQtn==1));

to_plot{4}=abs(rad_input.deltaZbaseline(rad_input.isQtn==1));

subplot(2,8,15)
easy_box(to_plot)
ylim([-0.1 1])
temp_labels={'buffered','potentiated','found in no rad','found in rad'};
xticks(1:length(temp_labels))
xticklabels(temp_labels)
xtickangle(45)
title('effect size in no rad')
ylabel('\DeltaZ')



clear to_plot
to_plot{1}=abs(hsp90_input.deltaZbuffer(logical(hsp90_qtn_idx.*v_buffer)));
to_plot{2}=abs(hsp90_input.deltaZbuffer(logical(hsp90_qtn_idx.*v_potentiate)));

to_plot{3}=abs(no_rad_input.deltaZbuffer(no_rad_input.isQtn==1));

to_plot{4}=abs(rad_input.deltaZbuffer(rad_input.isQtn==1));

subplot(2,8,16)
easy_box(to_plot)
ylim([-0.1 1])
temp_labels={'buffered','potentiated','found in no rad','found in rad'};
xticks(1:length(temp_labels))
xticklabels(temp_labels)
xtickangle(45)
title('effect size in rad')
ylabel('\DeltaZ')


end


