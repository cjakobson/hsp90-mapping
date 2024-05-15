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

%filebase='/Users/cjakobson/';
filebase='/Users/christopherjakobson/';

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
plot_rictor_reconstruction([3 4],dependency_directory,output_directory)


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



%Figure S1
figure('units','normalized','outerposition',[0 0 1 1])

%A
%v-Src assay
%from plotVsrc.m
subplot(2,8,1)
plot_vsrc(dependency_directory,output_directory)


%B
%scatter min glc -/+rad
subplot(2,4,2)
plot_trait_scatter('72h min glc_2%-rad','72h min glc_2%-rad',dependency_directory,output_directory)


%C
%n/a


%D
%heat shock mRNAs
%from analyzeRnaSeq.m
plot_heat_shock_rnas(dependency_directory,output_directory)


%E
%cross plot traits
subplot(2,4,5)
plot_trait_scatter('72h fluconazole_100uM_delta','72h rapamycin_5uM_delta',dependency_directory,output_directory)


%F
%n/a


%G
%effect size -- buffered vs potentiated
subplot(2,8,11)
plot_effect_size_interaction(dependency_directory,output_directory)


%H
%effect size -- buffered vs potentiated split by inverted/not
subplot(2,8,12)
plot_effect_size_line_crossing(dependency_directory,output_directory)



%I
%effect size histogam -- no rad vs Hsp90-dependent
subplot(2,4,7)
plot_effect_size_histogram(dependency_directory,output_directory)




%J
%Hsp90 vs linear effect size in no rad and rad
plot_effect_size_across_mapping(dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S1_1'],'-dsvg','-r0')
print([output_directory 'figure_S1_1'],'-djpeg','-r300')



figure('units','normalized','outerposition',[0 0 1 1])

%K
%RICTOR reconstruction in glucose
plot_rictor_reconstruction([1 2],dependency_directory,output_directory)


%L
%cumulative plot for effect size
subplot(2,4,2)
plot_cumulative_effect_size(dependency_directory,output_directory)


%M
%effect size by variant type
subplot(2,4,3)
plot_effect_size_type(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S1_2'],'-dsvg','-r0')
print([output_directory 'figure_S1_2'],'-djpeg','-r300')




toc

