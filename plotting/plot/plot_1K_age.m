function [] = plot_1K_age(plot_offset,dependency_directory,output_directory)



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


%first plot n_muts and AF as f(gene_age)
gene_age_info=readtable([dependency_directory 'Supplementary_Data_4_Doughty_et_al_2020.xlsx']);

group_labels={'Group I','Group II','Group III','WGD','Group IV','Group V'};

gene_age=nan(length(all_genes),1);
for i=1:length(all_genes)

    temp_idx=ismember(gene_age_info.Var1,all_genes{i});

    if sum(temp_idx)>0

        gene_age(i)=find(ismember(group_labels,gene_age_info.Var2{temp_idx}));

    end

end

clear to_plot
m=1;
for i=1:length(group_labels)
    
    temp_idx=gene_age==i;

    to_plot{m}=mis_n_muts(temp_idx);
    m=m+1;

end

subplot(2,8,plot_offset+1)
hold on
easy_box(to_plot)
ylim([0 0.2])
xlim([0.5 length(group_labels)+0.5])
xticks(1:length(group_labels))
xtickangle(45)
xticklabels(group_labels)
ylabel('missense variants per bp')
%axis square
for i=1:(length(to_plot)-1)
    [p h]=ranksum(to_plot{i},to_plot{end});
    text(i,0.18,num2str(p))
end


clear to_plot
m=1;
for i=1:length(group_labels)
    
    temp_idx=gene_age==i;

    to_plot{m}=reg_n_muts(temp_idx);
    m=m+1;

end

subplot(2,8,plot_offset+2)
hold on
easy_box(to_plot)
ylim([0 0.3])
xlim([0.5 length(group_labels)+0.5])
xticks(1:length(group_labels))
xtickangle(45)
xticklabels(group_labels)
ylabel('regulatory variants per bp')
%axis square
for i=1:(length(to_plot)-1)
    [p h]=ranksum(to_plot{i},to_plot{end});
    text(i,0.18,num2str(p))
end


end



