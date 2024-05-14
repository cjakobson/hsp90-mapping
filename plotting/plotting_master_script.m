%produce various plots for Hsp90 manuscript

clear

tic

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

filebase='/Users/cjakobson/';
%filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/hsp90-mapping/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/hsp90mapping/hsp90-mapping-dependencies/'];
output_directory=[filebase 'Dropbox/JaroszLab/hsp90mapping/manuscript-plots/'];

addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/parse'])
addpath([code_directory 'plotting/calculate'])
addpath([code_directory 'plotting/plot'])
addpath([code_directory 'data-prep'])

%Figure 1
figure('units','normalized','outerposition',[0 0 1 1])

%A
%n/a

%B
%rapamycin example segregant data
%from Hsp90mappingPlotsDelta.m
subplot(2,2,1)
plot_segregants_ranked('72h rapamycin_5uM',dependency_directory,output_directory)


%C
%LOD plot for rapamycin
%from plotQtnScores.m
subplot(4,1,3)
plot_lod('72h rapamycin_5uM',dependency_directory,output_directory)




%D
%RICTOR allele effect
%from hsp90examplesDelta.m
subplot(2,4,3)
plot_allele_effect(3176,'72h rapamycin_5uM',dependency_directory,output_directory)


%E
%LOD and QTN score for RICTOR
%from plotQtnScores.m
subplot(2,4,4)
plot_lod_qtn(3176,'72h rapamycin_5uM',dependency_directory,output_directory)


%F
%n/a


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_1_1'],'-dsvg','-r0')
print([output_directory 'figure_1_1'],'-djpeg','-r300')


figure('units','normalized','outerposition',[0 0 1 1])

%G
%CRISPR reconstruction of RICTOR
%from plotRadPhenotpying5b.m
%subplot(2,4,1)
plot_rictor_reconstruction(dependency_directory,output_directory)


%H
%RICTOR allele effect in 1K
%from analyze1K.m
subplot(2,8,3)
plot_1K_effect(3176,'rapamycin','48h',dependency_directory,output_directory)



%I
%cumulative frequency by pVal and molecular type
subplot(2,4,3)
plot_cumulative_pval(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_1_2'],'-dsvg','-r0')
print([output_directory 'figure_1_2'],'-djpeg','-r300')






%couple of exploratory plots for DFJ
figure('units','normalized','outerposition',[0 0 1 1])

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

subplot(2,8,1)
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

subplot(2,8,2)
easy_box(to_plot)
ylim([-0.1 1])
temp_labels={'buffered','potentiated','found in no rad','found in rad'};
xticks(1:length(temp_labels))
xticklabels(temp_labels)
xtickangle(45)
title('effect size in rad')
ylabel('\DeltaZ')


%cumulative plot for effect size
subplot(2,4,2)
plot_cumulative_effect_size(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_1_exploratory'],'-dsvg','-r0')
print([output_directory 'figure_1_exploratory'],'-djpeg','-r300')



toc

