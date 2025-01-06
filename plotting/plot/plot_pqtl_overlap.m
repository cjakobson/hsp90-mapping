function [] = plot_no_stress_overlap(dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

%just QTNs
input_data=input_data(input_data.isQtn==1,:);

conditions=unique(input_data.condition);

pqtl_data=readtable([dependency_directory 'linearPqtlOd_FDR_0.1.csv']);

%just HSC82/HSP82
idx_to_use=logical(ismember(pqtl_data.protein,'YMR186W')+...
    ismember(pqtl_data.protein,'YPL240C'));

pqtl_loci=unique(pqtl_data.index(idx_to_use))


dist_thresh=5;

for i=1:length(conditions)
    
    temp_loci=unique(input_data.index(ismember(input_data.condition,conditions{i})));
    
    clear min_dist
    for j=1:length(temp_loci)
        
        min_dist(j)=min(abs(pqtl_loci-temp_loci(j)));
        
    end
    
    f_overlap(i)=sum(min_dist<dist_thresh)/length(min_dist);

    n_overlap(i)=sum(min_dist<dist_thresh);
    
end

bar(f_overlap)
xticks(1:length(conditions))
xtickangle(45)
xticklabels(conditions)
xlim([0 length(conditions)+1])
title('all loci')
axis square
ylabel('fraction overlapping with Hsp90 pQTLs')

sum(n_overlap)

end

