function []= plot_trait_scatter(condition_to_plot1,condition_to_plot2,dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

load([dependency_directory 'radFilename.mat'])
load([dependency_directory 'radTrait.mat'])


trait_idx1=find(ismember(filename,condition_to_plot1));
trait_idx2=find(ismember(filename,condition_to_plot2));

v1=trait{trait_idx1};
v2=trait{trait_idx2};

hold on
scatter(v1,v2,5,'k','filled')
xlim([-3 3])
ylim(xlim)
axis square
[r p]=corr(v1,v2,'rows','complete');
text(1,-2,num2str(r))


end


