function [] = plot_1K_telomere(plot_offset,dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

gene_info=readtable([dependency_directory 'TableS4.xls']);


af_input=readtable([dependency_directory '1K_frequencies.csv']);

mis_n_muts=af_input.mis_n_muts;
reg_n_muts=af_input.reg_n_muts;

all_genes=af_input.Var1;

v_dist = calculate_telomere_distance_gene(dependency_directory);


telo_dist=nan(length(all_genes),1);
for i=1:length(all_genes)

    temp_idx=ismember(gene_info.Name,all_genes{i});

    if sum(temp_idx)>0

        telo_dist(i)=v_dist(temp_idx);

    end

end



%do these by bins instead as in existing figure
telo_bins=[0 0.05 0.25 0.5];



subplot(2,8,plot_offset+1)
hold on
v1=telo_dist;
v2=mis_n_muts;
clear to_plot
m=1;
for i=1:(length(telo_bins)-1)
    temp_idx=logical((v1>telo_bins(i)).*(v1<=telo_bins(i+1)));
    
    to_plot{m}=v2(temp_idx);
    m=m+1;
end
easy_box(to_plot)
ylim([0 0.2])
xlim([0.5 length(to_plot)+0.5])
%axis square
ylabel('missense variants per bp')
xlabel('position relative to telomere')
for i=2:length(to_plot)
    [p h]=ranksum(to_plot{1},to_plot{i});
    text(i,0.18,num2str(p))
end


subplot(2,8,plot_offset+2)
hold on
v1=telo_dist;
v2=reg_n_muts;
clear to_plot
m=1;
for i=1:(length(telo_bins)-1)
    temp_idx=logical((v1>telo_bins(i)).*(v1<=telo_bins(i+1)));
    
    to_plot{m}=v2(temp_idx);
    m=m+1;
end
easy_box(to_plot)
ylim([0 0.3])
xlim([0.5 length(to_plot)+0.5])
%axis square
ylabel('regulatory variants per bp')
xlabel('position relative to telomere')
for i=2:length(to_plot)
    [p h]=ranksum(to_plot{1},to_plot{i});
    text(i,0.18,num2str(p))
end



end



