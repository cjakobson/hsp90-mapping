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
%addpath([code_directory 'plotting/calculate'])
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
%from 




