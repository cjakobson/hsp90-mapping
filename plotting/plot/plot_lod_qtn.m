function [] = plot_lod_qtn(locus_to_plot,condition_to_plot,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([dependency_directory 'radFilename.mat'])
load([dependency_directory 'radTrait.mat'])


trait_to_use=trait{ismember(filename,[condition_to_plot '_delta'])};

[lod_1D,perm_thresh] = calculate_lod(trait_to_use,10,dependency_directory,output_directory);


%also QTN score
mapping_data_to_load=[condition_to_plot '_delta.mat'];
load([dependency_directory 'linear-hsp90/' mapping_data_to_load])

%assume input locus doesn't hav 
idx_to_use=find(ismember(posToMap,locus_to_plot+41));

if length(idx_to_use)>0
    temp_mat=ph2{idx_to_use};
    hold on
    temp_mat(temp_mat==-1)=NaN;
    v_to_plot=min(real(-log10(temp_mat)),[],2);

    v_to_plot(isnan(v_to_plot))=0;

    yyaxis right
    scatter(1:length(v_to_plot),v_to_plot,25,'k','filled')
    plot(v_to_plot,'-k','LineWidth',1)
    title([condition_to_plot ' ' num2str(locus_to_plot)])
    xlim([0 21])
    plot([11 11],ylim,':r')
    axis square
    ylabel('QTN score')
    
    yyaxis left
    left_lim=locus_to_plot-10;
    right_lim=locus_to_plot+10;
    
    lod_to_plot=lod_1D(left_lim:right_lim);
    
    scatter(1:length(lod_to_plot),lod_to_plot,25,'b','filled')
    plot(lod_to_plot,'--b','LineWidth',1)
    ylabel('LOD')
    
    

end




    
    
end
    