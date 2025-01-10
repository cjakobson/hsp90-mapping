function [] = plot_no_stress_overlap_sliding_p(dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


v_p_thresh=0:5:40;

input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

%just QTNs
input_data=input_data(input_data.isQtn==1,:);

for n=1:length(v_p_thresh)

    
    %filter on p
    temp_data=input_data(input_data.pVal>v_p_thresh(n),:);


    loci_min_glc=unique(temp_data.index(ismember(temp_data.condition,'min glc_2%')));

    conditions=unique(temp_data.condition);
    %remove min glc
    conditions(7)=[];

    dist_thresh=5;

    for i=1:length(conditions)

        temp_loci=unique(temp_data.index(ismember(temp_data.condition,conditions{i})));

        clear min_dist
        for j=1:length(temp_loci)

            min_dist(j)=min(abs(loci_min_glc-temp_loci(j)));

        end

        f_overlap(i)=sum(min_dist<dist_thresh)/length(min_dist);

    end

    mean_overlap(n)=mean(f_overlap);
    std_overlap(n)=std(f_overlap);
    
end

plot(mean_overlap,'k-')
%errorbar(mean_overlap,std_overlap,'k-')
xticks(1:length(v_p_thresh))
xtickangle(45)
xticklabels(v_p_thresh)
xlabel('p value threshold (-log_{10})')
ylabel('mean fraction overlapping')
axis square
ylim([0 1])



end

