function [] = plot_segregants_ranked(condition_to_plot,dependency_directory,output_directory)


    set(0,'DefaultLineLineWidth',1)
    set(0,'DefaultFigureColor','w')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultAxesLineWidth',1)
    
    blue=[43 172 226]./256;
    orange=[248 149 33]./256;
    grey=[128 128 128]./256;
    
    load([dependency_directory 'radFilename.mat'])
    load([dependency_directory 'radTrait.mat'])
    
    %colony size sorted histograms
    v1=trait{ismember(filename,[condition_to_plot '-rad'])};
    v2=trait{ismember(filename,[condition_to_plot '+rad'])};

    [v1_sort,sort_idx]=sort(v1,'ascend');
    v2_sort=v2(sort_idx);

    %subset
    v1_sort=v1_sort(1:100:length(v1_sort));
    v2_sort=v2_sort(1:100:length(v2_sort));

    hold on
    bar(v1_sort)
    bar(v2_sort)
    xlim([0 sum(~isnan(v1_sort))+1])
    ylim([-3 3])
    xlabel('segregants')
    ylabel('norm. growth')
    
    
end
    