function [] = plot_chaperone_enrichments(dependency_directory,output_directory)

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

hold on
clear v1 v2
for i=1:length(chaperone_query)
    temp_idx=ismember(all_genes,chaperone_query{i});
    v1(i)=relative_mat(1,temp_idx);
    v2(i)=interactor_pval(temp_idx);
end
bar(v1)
for i=1:length(v2)
    text(i,1.5+0.1*i,num2str(v2(i)))
end
xticks(1:length(chaperone_names))
xtickangle(45)
xticklabels(chaperone_names)
set(gca,'YScale','log')
ylim([0.5 3])
plot(xlim,[1 1],':k')
ylabel('rel. enrichment')



end