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
subplot(2,3,1)
plot_no_stress_overlap(dependency_directory,output_directory)


%extent of overlap with sliding p threshold
subplot(2,3,2)
plot_no_stress_overlap_sliding_p(dependency_directory,output_directory)



%overlap between Hsp90 QTNs and Hsc82/Hsp82 pQTLs
subplot(2,3,3)
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
        
        subplot(2,3,m)
        plot_trait_scatter(traits_to_plot{i},traits_to_plot{j},...
            dependency_directory,output_directory)
        m=m+1;
    
    end
    
end



set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R2_1'],'-dsvg','-r0')
print([output_directory 'figure_R2_1'],'-djpeg','-r300')




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
print([output_directory 'figure_R3_1'],'-dsvg','-r0')
print([output_directory 'figure_R3_1'],'-djpeg','-r300')




%some basic stats/analysis about synonymous

figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,3,1)
plot_variant_types(dependency_directory,output_directory)

%calculate total contribution by type

subplot(2,6,3)
plot_effect_size_ase_syn(dependency_directory,output_directory)


subplot(2,6,4)
plot_syn_delta_nte_ase(dependency_directory,output_directory)


subplot(2,3,3)
plot_syn_delta_nte_ase_hist(dependency_directory,output_directory)


subplot(2,6,7)
plot_syn_domain_boundaries(dependency_directory,output_directory)


% subplot(2,3,3)
% plot_syn_effect_size_ase(dependency_directory,output_directory)
% axis square


subplot(2,3,5)
plot_syn_pos_in_orf_ase(dependency_directory,output_directory)

subplot(2,3,6)
plot_syn_pos_in_orf_ase_nte(dependency_directory,output_directory)

%for FoldX -- need to include sign -- which allele is ref?




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R4_1'],'-dsvg','-r0')
print([output_directory 'figure_R4_1'],'-djpeg','-r300')


%various abundance correlations
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,3,1)
plot_effect_size_abundance(dependency_directory,output_directory)


%physical interactors
subplot(2,3,2)
plot_interactors_abundance(1,dependency_directory,output_directory)
title('physical')


subplot(2,3,3)
plot_interactors_abundance(2,dependency_directory,output_directory)
title('genetic')




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R5_1'],'-dsvg','-r0')
print([output_directory 'figure_R5_1'],'-djpeg','-r300')



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
print([output_directory 'figure_R6_1'],'-dsvg','-r0')
print([output_directory 'figure_R6_1'],'-djpeg','-r300')







%renormalize H2
figure('units','normalized','outerposition',[0 0 1 1])


subplot(2,3,1)
plot_heritability_explained(dependency_directory,output_directory)


subplot(2,3,2)
plot_heritability_explained_norm(dependency_directory,output_directory)


%fraction of modified vs all SNPs that are 90 vs 70 clients
subplot(2,3,3)
plot_gong_clients(dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R7_1'],'-dsvg','-r0')
print([output_directory 'figure_R7_1'],'-djpeg','-r300')




figure('units','normalized','outerposition',[0 0 1 1])


%plot Hsc82 and Hsp82 abundance vs growth in glc-/+rad
plot_hsp90_growth(0,dependency_directory,output_directory)




set(gcf,'PaperPositionMode','auto')
print([output_directory 'figure_R8_1'],'-dsvg','-r0')
print([output_directory 'figure_R8_1'],'-djpeg','-r300')





toc

