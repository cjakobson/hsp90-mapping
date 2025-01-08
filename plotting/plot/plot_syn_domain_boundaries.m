function [] = plot_syn_domain_boundaries(dependency_directory,output_directory)


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

[~,v_type]=variant_types(input_data.variantType);
v_type=v_type';

%only syn
input_data=input_data(v_type==2,:);


v_has_ase=logical(input_data.hasAseRad);

v_has_tag=logical(input_data.hasTagNoRad+input_data.hasTagRad);



%load domain boundaries
load([dependency_directory 'superfam_domain_boundaries.mat'])
load([dependency_directory 'gene_names.mat'])

min_boundary_dist=nan(height(input_data),1);
for i=1:height(input_data)

    temp_gene=input_data.gene1{i};
    temp_pos=ceil(str2num(input_data.encoded{i}(1:(end-3)))/3);

    gene_idx=find(ismember(genes_to_use,temp_gene));

    v_temp=domain_mat(gene_idx,:);
    boundary_idx=find(v_temp(2:end)~=v_temp(1:(end-1)));

    if ~isempty(boundary_idx)
    
        min_boundary_dist(i)=min(abs(boundary_idx-temp_pos));

    end

end


to_plot{1}=min_boundary_dist(logical(v_has_ase));
to_plot{2}=min_boundary_dist(logical(v_has_tag.*~v_has_ase));

%all other syn

background_data=readtable([dependency_directory 'variantInfoStructure.csv']);

[~,v_type_all]=variant_types(background_data.variantType);
v_type_all=v_type_all';

%only syn
background_data=background_data(v_type_all==2,:);

min_boundary_dist_all=nan(height(input_data),1);
for i=1:height(input_data)

    temp_gene=background_data.gene1{i};
    temp_pos=ceil(str2num(background_data.encoded{i}(1:(end-3)))/3);

    gene_idx=find(ismember(genes_to_use,temp_gene));

    v_temp=domain_mat(gene_idx,:);
    boundary_idx=find(v_temp(2:end)~=v_temp(1:(end-1)));

    if ~isempty(boundary_idx)
    
        min_boundary_dist_all(i)=min(abs(boundary_idx-temp_pos));

    end

end

to_plot{3}=min_boundary_dist_all;



easy_box(to_plot)
%ylim([0 0.5])
for i=1:length(to_plot)
    text(i,0,num2str(length(to_plot{i})))
end
for i=2:length(to_plot)
    [h p]=ttest2(to_plot{1},to_plot{i});
    text(i-0.5,600,num2str(p))
end
ylabel('distance to nearest domain boundary')
xticks(1:length(to_plot))
xtickangle(45)
xticklabels({'has ASE','no ASE','all segr. syn.'})
title('synonymous variants')
set(gca,'YScale','log')

end


