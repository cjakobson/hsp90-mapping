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
variantInfo=readtable([dependency_directory 'variantInfoStructure.csv']);
vDiff=variantInfo.chr(2:end)~=variantInfo.chr(1:(end-1));
breaks_to_plot=find(vDiff);


model_genotypes = parse_genotypes(dependency_directory,output_directory);
    

trait_to_use=trait{ismember(filename,[condition_to_plot '_delta'])};

%coarse LOD scoring to start (no adjustments) to ID hotspots
lod_1D=nan(1,length(model_genotypes(1,:)));
for q = 1:length(model_genotypes(1,:))
    r = corr(model_genotypes(:,q),trait_to_use,'rows','complete');
    lod_1D(q) = -length(trait_to_use)*log(1-r^2)/(2*log(10));
end

hold on
scatter(1:length(lod_1D),lod_1D,10,'k','filled')
xlim([1 length(lod_1D)])
ylim([0 200])
for k=1:length(breaks_to_plot)
    plot([breaks_to_plot(k) breaks_to_plot(k)],ylim,':k','LineWidth',1)
end




    
    
end
    