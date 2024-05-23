function [] = plot_lod(condition_to_plot,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

load([dependency_directory 'radFilename.mat'])
load([dependency_directory 'radTrait.mat'])


%to annotate chrs on plots
variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);
v_diff=variant_info.chr(2:end)~=variant_info.chr(1:(end-1));
breaks_to_plot=find(v_diff);



trait_to_use=trait{ismember(filename,[condition_to_plot '_delta'])};

[lod_1D,perm_thresh] = calculate_lod(trait_to_use,10,dependency_directory,output_directory);

hold on
scatter(1:length(lod_1D),lod_1D,10,'k','filled')
xlim([1 length(lod_1D)])
ylim([0 200])
for k=1:length(breaks_to_plot)
    plot([breaks_to_plot(k) breaks_to_plot(k)],ylim,':k','LineWidth',1)
end
plot(xlim,[perm_thresh perm_thresh],':r')
xlabel('genomic coordinate')
ylabel('LOD')



    
    
end
    