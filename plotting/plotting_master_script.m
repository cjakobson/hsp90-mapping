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
plot_trait_scatter('72h min glc_2%-rad','72h min glc_2%+rad',dependency_directory,output_directory)


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


close all


%Figure 2
figure('units','normalized','outerposition',[0 0 1 1])


%A
%ERG11 QTN vs LOD
subplot(2,4,1)
plot_lod_qtn(4975,'72h fluconazole_100uM',dependency_directory,output_directory)


%B
%ERG11 reconstruction in flc
%from readSgaDataErg11rad.m
plot_erg11_reconstruction(1,[1 2],2,dependency_directory,output_directory)


%C
%ASE -/+rad
subplot(2,4,3)
plot_ase(dependency_directory,output_directory)


%D
%effect size as function of distance from TSS
%from hsp90regulatoryPlots.m
subplot(2,4,4)
plot_tss_effect_size(dependency_directory,output_directory)



%E
%ERG11 ASE
subplot(2,8,9)
plot_tag_counts(4972,dependency_directory,output_directory)


subplot(2,8,10)
plot_allele_ratio(4972,dependency_directory,output_directory)



%F
%as D but for TES
subplot(2,4,6)
plot_tes_effect_size(dependency_directory,output_directory)



%G
%NFS1 reconstruction
%from readSgaDataNfs1rad.m
plot_nfs1_reconstruction([3 8],24,dependency_directory,output_directory)

plot_nfs1_reconstruction([2 7],26,dependency_directory,output_directory)


%H
%effect size of reg variants by interaction type and -/+ASE
subplot(2,4,8)
plot_reg_effect_size_ase(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2'],'-dsvg','-r0')
print([output_directory 'figure_2'],'-djpeg','-r300')



%Figure S2
figure('units','normalized','outerposition',[0 0 1 1])


%A
%ERG11 reconstruction in teb
%from readSgaDataErg11rad.m
plot_erg11_reconstruction(1,[3 4],0,dependency_directory,output_directory)


%B
%ERG11 reconstruction in glc
%from plotRadPhenotyping2.m
plot_erg11_reconstruction_glucose(1,2,dependency_directory,output_directory)


%C
%n/a


%D
%ASE power calculation
%from calculateASE.m
subplot(2,4,3)
plot_ase_power(dependency_directory,output_directory)



%E
%ASE -/+rad for all sites
subplot(2,4,4)
plot_ase_all(dependency_directory,output_directory)



%F
%effect size for has ASE/not for all SNPs
subplot(2,8,9)
plot_effect_size_ase(dependency_directory,output_directory)



%G
%effect size of protein-coding variants by interaction type and -/+ASE
subplot(2,4,6)
plot_coding_effect_size_ase(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S2'],'-dsvg','-r0')
print([output_directory 'figure_S2'],'-djpeg','-r300')



close all



%Figure 3
figure('units','normalized','outerposition',[0 0 1 1])
%A
%n/a


%B
%chaperone and cochaperone enrichments
subplot(2,4,1)
plot_chaperone_enrichments(dependency_directory,output_directory)


%C and C inset -- zoom in
%fraction of modified loci interacting
plot_fraction_interactors(1,dependency_directory,output_directory)



%D
%scatter plot fraction enriched/enrichment pVal
subplot(2,4,4)
plot_enrichment_scatter(dependency_directory,output_directory)



%E
%kinase and TF enrichments
%from kinomeAnalysis.m
subplot(2,2,3)
plot_kinase_tf_enrichments(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_3'],'-dsvg','-r0')
print([output_directory 'figure_3'],'-djpeg','-r300')



%Figure S3
figure('units','normalized','outerposition',[0 0 1 1])

%A
%enrichment for all kinases
subplot(2,2,1)
plot_all_kinases(dependency_directory,output_directory)

%B
%also for all other TFs?
subplot(2,2,2)
plot_all_tfs(dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S3'],'-dsvg','-r0')
print([output_directory 'figure_S3'],'-djpeg','-r300')



close all



%Figure 4
figure('units','normalized','outerposition',[0 0 1 1])

%A
%n/a

%B
%n/a


%C 
%mutational step effect sizes
loci_to_plot(1)=4835;    %IMA1/MAL1
loci_to_plot(2)=886;     %MAL3
loci_to_plot(3)=426;     %GAL10

ref_allele=[1 1 -1];

subplot(2,8,1)
plot_mutational_steps(loci_to_plot,ref_allele,'72h min mal_2%',...
    dependency_directory,output_directory)


subplot(2,8,2)
plot_mutational_steps(loci_to_plot,ref_allele,'72h min raf_2%',...
    dependency_directory,output_directory)


%D
%same for flc/glc
loci_to_plot(1)=4974;    %ERG11
loci_to_plot(2)=4292;     %PDR1
loci_to_plot(3)=10060;     %GRE2

ref_allele=[1 -1 1];

subplot(2,8,3)
plot_mutational_steps(loci_to_plot,ref_allele,'72h fluconazole_100uM',...
    dependency_directory,output_directory)


subplot(2,8,4)
plot_mutational_steps(loci_to_plot,ref_allele,'72h min glc_2%',...
    dependency_directory,output_directory)


%E
%n/a


%F
%gene age for interacting variants
subplot(2,4,3)
plot_gene_age_all(dependency_directory,output_directory)

%G
%only regulatory
subplot(2,4,4)
plot_gene_age_type(3,dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_4'],'-dsvg','-r0')
print([output_directory 'figure_4'],'-djpeg','-r300')





%Figure S4
figure('units','normalized','outerposition',[0 0 1 1])

%A
%heritability explained
subplot(2,4,1)
plot_heritability_explained(dependency_directory,output_directory)


%B
%modified QTLs per trait
subplot(2,4,2)
plot_qtls_per_trait(dependency_directory,output_directory)


%C
%maltose allele effect
subplot(2,4,3)
plot_allele_effect(4835,'72h min mal_2%',dependency_directory,output_directory)




%D
loci_to_plot(1)=4835;    %IMA1/MAL1
loci_to_plot(2)=886;     %MAL3
loci_to_plot(3)=426;     %GAL10

ref_allele=[1 1 -1];

subplot(2,2,3)
plot_allele_combintations(loci_to_plot,ref_allele,'72h min mal_2%',...
    dependency_directory,output_directory)


%E
subplot(2,2,4)
plot_allele_combintations(loci_to_plot,ref_allele,'72h min raf_2%',...
    dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S4_1'],'-dsvg','-r0')
print([output_directory 'figure_S4_1'],'-djpeg','-r300')




figure('units','normalized','outerposition',[0 0 1 1])


%F
%modified QTN pleiotropy
subplot(2,4,1)
plot_allele_pleiotropy_scatter(dependency_directory,output_directory)




%new supp panel -- missense gene age
%G
subplot(2,4,2)
plot_gene_age_type(1,dependency_directory,output_directory)


%H
%gene age for buffered
subplot(2,4,3)
plot_gene_age_interaction(1,dependency_directory,output_directory)



%new panel
%also for potentiated
%I
subplot(2,4,4)
plot_gene_age_interaction(2,dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S4_2'],'-dsvg','-r0')
print([output_directory 'figure_S4_2'],'-djpeg','-r300')


close all


%Figure 5
figure('units','normalized','outerposition',[0 0 1 1])

%A
%n/a


%B
%n/a


%C
%ERG11 CRISPRi experiment
%raw data for YJM975 flc and flc+rad
plot_erg11_crispri(0,dependency_directory,output_directory)


%D
%normalized YJM975 data in flc
subplot(2,4,3)
plot_erg11_crispri_norm(1,1,dependency_directory,output_directory)


%E
%effect size by haploinsuff, ess, etc
subplot(2,8,7)
plot_haplo_interaction(1,dependency_directory,output_directory)

subplot(2,8,8)
plot_haplo_interaction(2,dependency_directory,output_directory)


%F
%essential and haploinsufficient analyses for human
subplot(2,8,9)
plot_haplo_human(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_5'],'-dsvg','-r0')
print([output_directory 'figure_5'],'-djpeg','-r300')



%Figure S5
figure('units','normalized','outerposition',[0 0 1 1])

%A
%n/a


%B
%Erg11 variants in flc/rad
subplot(2,4,1)
plot_erg11_mutants(1,dependency_directory,output_directory)

%C
%same in teb
subplot(2,4,2)
plot_erg11_mutants(2,dependency_directory,output_directory)


%D
%YJM975 crispri in glucose (no drug)
subplot(2,4,3)
plot_erg11_crispri_norm_glc(1,1,dependency_directory,output_directory)


%E
%norm RM11 in flc
subplot(2,4,4)
plot_erg11_crispri_norm(2,1,dependency_directory,output_directory)


%F
%n/a


%G
subplot(2,8,9)
plot_haplo_all(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S5'],'-dsvg','-r0')
print([output_directory 'figure_S5'],'-djpeg','-r300')


close all


%Figure S6
figure('units','normalized','outerposition',[0 0 1 1])

plot_ase_reproducibility(dependency_directory,output_directory)

set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S6'],'-dsvg','-r0')
print([output_directory 'figure_S6'],'-djpeg','-r300')




%Figure S7
figure('units','normalized','outerposition',[0 0 1 1])

plot_tpm_reproducibility(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S7'],'-dsvg','-r0')
print([output_directory 'figure_S7'],'-djpeg','-r300')



close all


toc




