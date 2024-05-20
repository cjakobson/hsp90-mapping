function [] = plot_fraction_interactors(plot_offset,dependency_directory,output_directory)

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


load([output_directory 'biogrid_data.mat'])


mapping_input=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);
%only QTNs
mapping_input=mapping_input(mapping_input.isQtn==1,:);

input_genes{1}=mapping_input.gene1;
input_genes{1}(cellfun(@isempty,input_genes{1}))={'NA'};
input_genes{2}=mapping_input.gene2;
input_genes{2}(cellfun(@isempty,input_genes{2}))={'NA'};


variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

input_genes{3}=[variant_info.gene1;variant_info.gene2];
input_genes{3}(cellfun(@isempty,input_genes{3}))=[];

[overlap_mat,interactor_pval,relative_mat] = ...
    calculate_fraction_interactors(input_genes,dependency_directory,output_directory);

%Hsp90/70 machinery
%see https://journals.asm.org/doi/10.1128/mmbr.05018-11
chaperone_query={'YMR186W','YPL240C',...
    'YOR027W','YNL064C','YDR168W',...
    'YAL005C','YLL024C','YBL075C','YER103W',...
    'YDL229W','YNL209W','YPL106C','YBR169C','YGR123C',...
    'YKL117W','YJR032W','YDR214W','YOR057W'};
chaperone_names={'Hsc82','Hsp82',...
    'Sti1','Ydj1','Cdc37',...
    'Ssa1','Ssa2','Ssa3','Ssa4',...
    'Ssb1','Ssb2','Sse1','Sse2','Ppt1'...
    'Sba1','Cpr7','Aha1','Sgt1',};

subplot(2,4,plot_offset+1)
hold on
[v_sort,sort_idx]=sort(overlap_mat(1,:));
v_sort(v_sort<=1e-2)=1e-2;

scatter(1:length(v_sort),v_sort,10,'k','filled')
%label key players
for i=1:length(chaperone_query)
    temp_idx1=find(ismember(all_genes,chaperone_query{i}));
    temp_idx2=find(ismember(sort_idx,temp_idx1));
    scatter(temp_idx2,v_sort(temp_idx2),25,'r','filled')
    text(temp_idx2,v_sort(temp_idx2),chaperone_names{i})
end
xlim([0 length(sort_idx)])
axis square
set(gca,'YScale','log')
ylim([1e-2 1])


subplot(2,4,plot_offset+2)
n_to_output=100;

%inset
hold on
scatter(1:length(v_sort),v_sort,10,grey,'filled')
%label key players
for i=1:length(chaperone_query)
    temp_idx1=find(ismember(all_genes,chaperone_query{i}));
    temp_idx2=find(ismember(sort_idx,temp_idx1));
end
%how many of the others are themselves Hsp90 interactors?
hsp90_idx=logical(interaction_mat(ismember(all_labels,'HSC82'),:)+...
    interaction_mat(ismember(all_labels,'HSC82'),:));
hsp90_interactors=unique(all_genes(hsp90_idx));
m=0;
for i=(length(v_sort)-n_to_output+1):length(v_sort)
    if v_sort(i)>0.25
        text(i,v_sort(i),all_labels{sort_idx(i)})
    end
    if sum(ismember(hsp90_interactors,all_genes{sort_idx(i)}))>0
        scatter(i,v_sort(i),25,'k','filled')
        m=m+1;
    end
end
text(5800,0.9,[num2str(m) ' interact w Hsp90'])
xlim([length(sort_idx)-n_to_output+1 length(sort_idx)])
axis square
set(gca,'YScale','log')
ylim([1e-1 1])





end
