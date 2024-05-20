function [] = plot_qtls_per_trait(dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

conditions=unique(input_data.condition);
times=unique(input_data.time);

m=1;
for i=1:length(conditions)
    
    condition_idx=ismember(input_data.condition,conditions{i});
    
    for j=1:length(times)
        
        time_idx=ismember(input_data.time,times{j});
        
        n_qtls(m)=sum(condition_idx.*time_idx);
        m=m+1;
        
    end
    
end

hold on
histogram(n_qtls,[50:25:250])
axis square


end


