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
%output_directory=[filebase 'Dropbox/JaroszLab/hsp90mapping/manuscript-plots/'];
output_directory=[filebase 'Dropbox/JaroszLab/hsp90mapping/revision-plots/'];

addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/parse'])
addpath([code_directory 'plotting/calculate'])
addpath([code_directory 'plotting/plot'])
addpath([code_directory 'data-prep'])




%Figure R1
figure('units','normalized','outerposition',[0 0 1 1])

%overlap of min glc with other conditions
subplot(2,4,1)
plot_no_stress_overlap(dependency_directory,output_directory)


%extent of overlap with sliding p threshold
subplot(2,4,2)
plot_no_stress_overlap_sliding_p(dependency_directory,output_directory)



%overlap between Hsp90 QTNs and Hsc82/Hsp82 pQTLs
subplot(2,4,3)
plot_pqtl_overlap(dependency_directory,output_directory)






set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R1_1'],'-dsvg','-r0')
print([output_directory 'figure_R1_1'],'-djpeg','-r300')






%Figure R2
figure('units','normalized','outerposition',[0 0 1 1])

traits_to_plot={'72h min glc_2%-rad','72h min glc_2%+rad',...
    '72h rapamycin_5uM-rad','72h rapamycin_5uM+rad'};
    
m=1;
for i=1:length(traits_to_plot)
    
    for j=(i+1):length(traits_to_plot)
        
        subplot(2,4,m)
        plot_trait_scatter(traits_to_plot{i},traits_to_plot{j},...
            dependency_directory,output_directory)
        m=m+1;
    
    end
    
end



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R2_1'],'-dsvg','-r0')
print([output_directory 'figure_R2_1'],'-djpeg','-r300')


%renormalize H2
figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,4,1)
plot_heritability_explained(dependency_directory,output_directory)



plot_heritability_explained_norm(1,dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R3_1'],'-dsvg','-r0')
print([output_directory 'figure_R3_1'],'-djpeg','-r300')





%some basic stats/analysis about synonymous

figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,4,1)
plot_variant_types(dependency_directory,output_directory)

%calculate total contribution by type

subplot(2,8,5)
plot_effect_size_ase_syn(dependency_directory,output_directory)



subplot(2,4,4)
plot_syn_effect_size_ase(dependency_directory,output_directory)

% subplot(2,8,4)
% plot_syn_delta_cai_ase(dependency_directory,output_directory)

subplot(2,8,9)
plot_syn_delta_nte_ase(dependency_directory,output_directory)



% subplot(2,4,3)
% plot_syn_pos_in_orf_ase(dependency_directory,output_directory)






set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R4_1'],'-dsvg','-r0')
print([output_directory 'figure_R4_1'],'-djpeg','-r300')





%subset biogrid analysis

figure('units','normalized','outerposition',[0 0 1 1])

%physical only
subplot(2,4,1)
plot_chaperone_enrichments(1,1,dependency_directory,output_directory)
title('physical only')

% subplot(2,4,2)
% plot_enrichment_scatter(1,1,dependency_directory,output_directory)
% title('physical only')

%genetic only
subplot(2,4,2)
plot_chaperone_enrichments(2,1,dependency_directory,output_directory)
title('genetic only')

% subplot(2,4,4)
% plot_enrichment_scatter(2,1,dependency_directory,output_directory)
% title('genetic only')


%coding vs non-coding
subplot(2,4,3)
plot_chaperone_enrichments_type(0,1,1,dependency_directory,output_directory)
title('coding')


subplot(2,4,4)
plot_chaperone_enrichments_type(0,1,3,dependency_directory,output_directory)
title('regulatory')



%all four combos
subplot(2,4,5)
plot_chaperone_enrichments_type(1,1,1,dependency_directory,output_directory)
title('physical interactions -- coding')

subplot(2,4,6)
plot_chaperone_enrichments_type(2,1,1,dependency_directory,output_directory)
title('genetic interactions -- coding')


subplot(2,4,7)
plot_chaperone_enrichments_type(1,1,3,dependency_directory,output_directory)
title('physical interactions -- regulatory')

subplot(2,4,8)
plot_chaperone_enrichments_type(2,1,3,dependency_directory,output_directory)
title('genetic interactions -- regulatory')






set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R5_1'],'-dsvg','-r0')
print([output_directory 'figure_R5_1'],'-djpeg','-r300')



%various abundance correlations
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,4,1)
plot_effect_size_abundance(dependency_directory,output_directory)


%physical interactors
subplot(2,4,2)
plot_interactors_abundance(1,dependency_directory,output_directory)
title('physical')


subplot(2,4,3)
plot_interactors_abundance(2,dependency_directory,output_directory)
title('genetic')




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R6_1'],'-dsvg','-r0')
print([output_directory 'figure_R6_1'],'-djpeg','-r300')



%correlate various hub expression with Hsc82/Hsp82 amongst F6 haploids
hubs_to_correlate={'STI1','YDJ1','HEK2','ACT1','BFR1','DHH1','SSA1','SSB1'};
systematic_to_correlate={'YOR027W','YNL064C','YBL032W','YFL039C',...
    'YOR198C','YDL160C','YAL005C','YDL229W'};

figure('units','normalized','outerposition',[0 0 1 1])

for i=1:length(hubs_to_correlate)
    
    subplot(2,4,i)
    plot_hub_correlation(hubs_to_correlate{i},systematic_to_correlate{i},...
        dependency_directory,output_directory)
    
end



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R7_1'],'-dsvg','-r0')
print([output_directory 'figure_R7_1'],'-djpeg','-r300')



ynstbaevr



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
subplot(2,8,5)
plot_allele_effect(3176,'72h rapamycin_5uM',dependency_directory,output_directory)


%E
%RICTOR allele effect in 1K
%from analyze1K.m
subplot(2,8,6)
plot_1K_effect(3176,'rapamycin','48h',dependency_directory,output_directory)





set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_1_1'],'-dsvg','-r0')
print([output_directory 'figure_1_1'],'-djpeg','-r300')


figure('units','normalized','outerposition',[0 0 1 1])


%F
%cumulative plot for effect size
subplot(2,4,1)
plot_cumulative_effect_size(dependency_directory,output_directory)


%G
%cumulative frequency by pVal and molecular type
subplot(2,4,2)
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
%heritability explained
subplot(2,4,6)
plot_heritability_explained(dependency_directory,output_directory)



%G
%n/a


%H
%effect size -- buffered vs potentiated
subplot(2,8,13)
plot_effect_size_interaction(dependency_directory,output_directory)


%I
%effect size -- buffered vs potentiated split by inverted/not
subplot(2,8,14)
plot_effect_size_line_crossing(dependency_directory,output_directory)



%J
%effect size histogam -- no rad vs Hsp90-dependent
subplot(2,4,8)
plot_effect_size_histogram(dependency_directory,output_directory)






set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S1_1'],'-dsvg','-r0')
print([output_directory 'figure_S1_1'],'-djpeg','-r300')



figure('units','normalized','outerposition',[0 0 1 1])




%K
%Hsp90 vs linear effect size in no rad and rad
plot_effect_size_across_mapping(0,dependency_directory,output_directory)


plot_af_across_mapping(2,dependency_directory,output_directory)




%L
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
%LOD and QTN score for RICTOR
%from plotQtnScores.m
subplot(2,4,1)
plot_lod_qtn(3176,'72h rapamycin_5uM',dependency_directory,output_directory)


%B
%n/a


%C
%CRISPR reconstruction of RICTOR
%from plotRadPhenotpying5b.m
%subplot(2,4,1)
plot_rictor_reconstruction([3 4],3,dependency_directory,output_directory)


%D
%ERG11 QTN vs LOD
subplot(2,4,3)
plot_lod_qtn(4975,'72h fluconazole_100uM',dependency_directory,output_directory)


%E
%ERG11 reconstruction in flc
%from readSgaDataErg11rad.m
plot_erg11_reconstruction(1,[1 2],6,dependency_directory,output_directory)


%F
%ERG11 ASE
subplot(2,8,9)
plot_total_counts(4972,dependency_directory,output_directory)


subplot(2,8,10)
plot_tag_counts(4972,dependency_directory,output_directory)


subplot(2,8,11)
plot_allele_ratio(4972,dependency_directory,output_directory)






set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_2_1'],'-dsvg','-r0')
print([output_directory 'figure_2_1'],'-djpeg','-r300')






%Figure S2
figure('units','normalized','outerposition',[0 0 1 1])


%A
%RICTOR reconstruction in glucose
plot_rictor_reconstruction([1 2],1,dependency_directory,output_directory)


%B
%ERG11 allele effect
%from hsp90examplesDelta.m
subplot(2,8,3)
plot_allele_effect(4974,'72h fluconazole_100uM',dependency_directory,output_directory)


%C
%ERG11 reconstruction in teb
%from readSgaDataErg11rad.m
plot_erg11_reconstruction(1,[3 4],4,dependency_directory,output_directory)


%D
%ERG11 reconstruction in glc
%from plotRadPhenotyping2.m
plot_erg11_reconstruction_glucose(1,6,dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S2_1'],'-dsvg','-r0')
print([output_directory 'figure_S2_1'],'-djpeg','-r300')



close all



%Figure 3
figure('units','normalized','outerposition',[0 0 1 1])


%A
%ASE -/+rad
subplot(2,4,1)
plot_ase(dependency_directory,output_directory)


%B
%effect size as function of distance from TSS
%from hsp90regulatoryPlots.m
subplot(2,4,2)
plot_tss_effect_size(dependency_directory,output_directory)


%C
%as D but for TES
subplot(2,4,3)
plot_tes_effect_size(dependency_directory,output_directory)



%D
%NFS1 reconstruction
%from readSgaDataNfs1rad.m
plot_nfs1_reconstruction([3 8],16,dependency_directory,output_directory)

plot_nfs1_reconstruction([2 7],18,dependency_directory,output_directory)


%E
%effect size of reg variants by interaction type and -/+ASE
subplot(2,4,6)
plot_reg_effect_size_ase(dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_3_1'],'-dsvg','-r0')
print([output_directory 'figure_3_1'],'-djpeg','-r300')



%Figure S3
figure('units','normalized','outerposition',[0 0 1 1])

%A
%n/a


%B
%ASE power calculation
%from calculateASE.m
subplot(2,4,1)
plot_ase_power(dependency_directory,output_directory)



%C
%ASE -/+rad for all sites
subplot(2,4,2)
plot_ase_all(dependency_directory,output_directory)



%D
%effect size for has ASE/not for all SNPs
subplot(2,8,5)
plot_effect_size_ase(dependency_directory,output_directory)



%E
%effect size of protein-coding variants by interaction type and -/+ASE
subplot(2,4,4)
plot_coding_effect_size_ase(dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S3_1'],'-dsvg','-r0')
print([output_directory 'figure_S3_1'],'-djpeg','-r300')


%Figure 4
figure('units','normalized','outerposition',[0 0 1 1])


%A
%n/a


%B
%chaperone and cochaperone enrichments
subplot(2,4,1)
plot_chaperone_enrichments(0,1,dependency_directory,output_directory)


%C and C inset -- zoom in
%fraction of modified loci interacting
plot_fraction_interactors(1,0,1,dependency_directory,output_directory)



%D
%scatter plot fraction enriched/enrichment pVal
subplot(2,4,4)
plot_enrichment_scatter(0,1,dependency_directory,output_directory)



%E
%kinase and TF enrichments
%from kinomeAnalysis.m
subplot(2,2,3)
plot_kinase_tf_enrichments(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_4_1'],'-dsvg','-r0')
print([output_directory 'figure_4_1'],'-djpeg','-r300')



%Figure S4
figure('units','normalized','outerposition',[0 0 1 1])

%A
%all Hsp70s
subplot(2,2,1)
plot_chaperone_enrichments(0,2,dependency_directory,output_directory)



%B
%all Hsp40s
subplot(2,2,2)
plot_chaperone_enrichments(0,3,dependency_directory,output_directory)



%C
%enrichment for all kinases
subplot(2,1,2)
plot_all_kinases(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S4_1'],'-dsvg','-r0')
print([output_directory 'figure_S4_1'],'-djpeg','-r300')




figure('units','normalized','outerposition',[0 0 1 1])

%D
%also for all other TFs
subplot(2,1,1)
plot_all_tfs(dependency_directory,output_directory)




%E
%Ub ligases
%see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3454868/
subplot(2,1,2)
plot_chaperone_enrichments(1,4,dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S4_2'],'-dsvg','-r0')
print([output_directory 'figure_S4_2'],'-djpeg','-r300')





close all



%Figure 5
figure('units','normalized','outerposition',[0 0 1 1])


%A
loci_to_plot(1)=4835;    %IMA1/MAL1
loci_to_plot(2)=886;     %MAL3
loci_to_plot(3)=426;     %GAL10

ref_allele=[1 1 -1];

subplot(2,2,1)
plot_allele_combintations(loci_to_plot,ref_allele,'72h min mal_2%',...
    dependency_directory,output_directory)


subplot(2,2,2)
plot_allele_combintations(loci_to_plot,ref_allele,'72h min raf_2%',...
    dependency_directory,output_directory)



%B
%n/a

%C
%mutational step effect sizes
loci_to_plot(1)=4835;    %IMA1/MAL1
loci_to_plot(2)=886;     %MAL3
loci_to_plot(3)=426;     %GAL10

ref_allele=[1 1 -1];

subplot(2,8,9)
condition_to_plot='72h min mal_2%';
plot_mutational_steps(loci_to_plot,ref_allele,condition_to_plot,...
    dependency_directory,output_directory)


subplot(2,8,10)
condition_to_plot='72h min raf_2%';
plot_mutational_steps(loci_to_plot,ref_allele,condition_to_plot,...
    dependency_directory,output_directory)

%D
%n/a

%E
%same for flc/glc
loci_to_plot(1)=4974;    %ERG11
loci_to_plot(2)=4292;     %PDR1
loci_to_plot(3)=10060;     %GRE2

ref_allele=[1 -1 1];

subplot(2,8,11)
condition_to_plot='72h fluconazole_100uM';
plot_mutational_steps(loci_to_plot,ref_allele,condition_to_plot,...
    dependency_directory,output_directory)


subplot(2,8,12)
condition_to_plot='72h min glc_2%';
plot_mutational_steps(loci_to_plot,ref_allele,condition_to_plot,...
    dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_5_1'],'-dsvg','-r0')
print([output_directory 'figure_5_1'],'-djpeg','-r300')



%Figure 6
figure('units','normalized','outerposition',[0 0 1 1])


%A
%n/a


%B
%gene age for interacting variants
subplot(2,4,1)
plot_gene_age_all(dependency_directory,output_directory)
title('all QTNs')

%C
%only regulatory
subplot(2,4,2)
plot_gene_age_type(3,dependency_directory,output_directory)
title('regulatory')


%D
%density of interactions vs genes
subplot(2,4,3)
plot_gene_density_telomere(dependency_directory,output_directory)


%E
%polymorphism by gene age
plot_1K_age(8,dependency_directory,output_directory)


%F
%polymorphism in telomere
plot_1K_telomere(10,dependency_directory,output_directory)



%G
%buffered vs potentiated effect size by telomere distance
subplot(2,4,4)
plot_effect_size_telomere(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_6_1'],'-dsvg','-r0')
print([output_directory 'figure_6_1'],'-djpeg','-r300')




%Figure S5
figure('units','normalized','outerposition',[0 0 1 1])

%A

%B
%modified QTLs per trait
subplot(2,4,1)
plot_qtls_per_trait(dependency_directory,output_directory)


%C
%maltose allele effect
subplot(2,4,2)
plot_allele_effect(4835,'72h min mal_2%',dependency_directory,output_directory)


%C
%modified QTN pleiotropy
subplot(2,4,3)
plot_allele_pleiotropy_scatter(dependency_directory,output_directory)


%D
%n/a


%new supp panel -- missense gene age
%E
subplot(2,4,5)
plot_gene_age_type(1,dependency_directory,output_directory)
title('missense')


%F
%gene age for buffered
subplot(2,4,6)
plot_gene_age_interaction(1,dependency_directory,output_directory)
title('all buffered')


%new panel
%also for potentiated
%G
subplot(2,4,7)
plot_gene_age_interaction(2,dependency_directory,output_directory)
title('all potentiated')



%H
%density of interactions vs variants
subplot(2,4,8)
plot_density_telomere_type(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S5_1'],'-dsvg','-r0')
print([output_directory 'figure_S5_1'],'-djpeg','-r300')




close all




%Figure 7
figure('units','normalized','outerposition',[0 0 1 1])


%A
%n/a


%B
%Erg11 variants in flc/rad
subplot(2,4,1)
plot_erg11_mutants(1,dependency_directory,output_directory)
title('fluconazole')



%C
%n/a


%D
%n/a


%E
%ERG11 CRISPRi experiment
%raw data for YJM975 flc and flc+rad
plot_erg11_crispri(1,dependency_directory,output_directory)



%F
%normalized YJM975 data in flc
subplot(2,4,5)
plot_erg11_crispri_norm(1,1,dependency_directory,output_directory)


%G
%effect size by haploinsuff, ess, etc
subplot(2,8,11)
plot_haplo_interaction(1,dependency_directory,output_directory)

subplot(2,8,12)
plot_haplo_interaction(2,dependency_directory,output_directory)


%H
%essential and haploinsufficient analyses for human
subplot(2,8,13)
plot_haplo_human(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_7_1'],'-dsvg','-r0')
print([output_directory 'figure_7_1'],'-djpeg','-r300')



%Figure S6
figure('units','normalized','outerposition',[0 0 1 1])

%A
%same in teb
subplot(2,4,1)
plot_erg11_mutants(2,dependency_directory,output_directory)
title('tebuconazole')


%B
%YJM975 crispri in glucose (no drug)
subplot(2,4,2)
plot_erg11_crispri_norm_glc(1,1,dependency_directory,output_directory)


%C
%norm RM11 in flc
subplot(2,4,3)
plot_erg11_crispri_norm(2,1,dependency_directory,output_directory)


%D
%n/a


%E
subplot(2,8,7)
plot_haplo_all(dependency_directory,output_directory)


set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S6_1'],'-dsvg','-r0')
print([output_directory 'figure_S6_1'],'-djpeg','-r300')


close all


%Figure S7
figure('units','normalized','outerposition',[0 0 1 1])

plot_ase_reproducibility(dependency_directory,output_directory)

set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S7_1'],'-dsvg','-r0')
print([output_directory 'figure_S7_1'],'-djpeg','-r300')


figure('units','normalized','outerposition',[0 0 1 1])

plot_tpm_reproducibility(dependency_directory,output_directory)



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_S7_2'],'-dsvg','-r0')
print([output_directory 'figure_S7_2'],'-djpeg','-r300')



close all


toc




