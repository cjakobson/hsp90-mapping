function [] = plot_haplo_human(dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


%load pLOF data
lof_table=readtable([dependency_directory 'supplementary_dataset_11_full_constraint_metrics.tsv'],...
    'FileType','text');

idx_to_keep=ismember(lof_table.canonical,'true');
lof_table=lof_table(idx_to_keep,:);

lof_table(isnan(lof_table.oe_lof_upper_rank),:)=[];

[~,sort_idx]=sort(lof_table.oe_lof_upper_rank,'descend');



%load haploinsuff genes
hap_table=readtable([dependency_directory '1-s2.0-S2211124722003175-mmc2.xlsx']);
hap_genes=hap_table.Symbol;



n_bins=5;

%most to least essential quintiles
quintile_size=floor(length(sort_idx)/n_bins);

for i=1:n_bins
    
    temp_idx=(i-1)*quintile_size+(1:quintile_size);
    gene_sets{i}=lof_table.gene(sort_idx(temp_idx));
    
end



client_table=readtable([dependency_directory 'picardHsp90Interactors.csv']);

to_keep=logical(ismember(client_table.Interactor_A,{'HSP90AA1','HSP90BB1'})+...
    ismember(client_table.Interactor_B,{'HSP90AA1','HSP90BB1'}));

client_table=client_table(to_keep,:);

%only keep the ones that arent hsp90
v_clients=unique([client_table.Interactor_A; client_table.Interactor_B]);
v_clients(ismember(v_clients,{'HSP90AA1','HSP90BB1'}))=[];


for i=1:length(gene_sets)
    
    n_int(i)=sum(ismember(gene_sets{i},v_clients));
    
    f_int(i)=n_int(i)/length(gene_sets{i});
    
end

n_int(i+1)=sum(ismember(hap_genes,v_clients));
f_int(i+1)=n_int(i+1)/length(hap_genes);


bar(f_int)

ylabel('fraction interactors')



end

