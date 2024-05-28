function [] = plot_kinase_tf_enrichments(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

qtn_idx=input_data.isQtn==1;
temp_genes=unique([input_data.gene1(qtn_idx);input_data.gene2(qtn_idx)]);
%temp_genes=[input_data.gene1(qtn_idx);input_data.gene2(qtn_idx)];
temp_genes(cellfun(@isempty,temp_genes))=[];

qtn_genes=temp_genes;


variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

temp_genes=unique([variant_info.gene1;variant_info.gene2]);
temp_genes(cellfun(@isempty,temp_genes))=[];

all_segregating_genes=temp_genes;
all_segregating_genes(ismember(all_segregating_genes,qtn_genes))=[];




v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);


%kinase data from MOTIPS/https://www.science.org/doi/10.1126/scisignal.2000482
%(dataset 2)

temp_names=dir([dependency_directory 'Dataset S2']);

m=1;
for i=4:length(temp_names)
    
    kinases{m}=temp_names(i).name;
    m=m+1;
    
end

kinases=kinases';

%threshold on likelihood to include
thresh=1;
for i=1:length(kinases)
    
    temp_names=dir([dependency_directory 'Dataset S2/' kinases{i}]);
    to_get=temp_names(3).name;
    
    temp_data=readtable([dependency_directory 'Dataset S2/' kinases{i} '/' to_get],'FileType','text');
    
    targets{i}=temp_data.ORF(temp_data.Likelihood>=thresh);
    
    target_overlap{i,1}=intersect(targets{i},qtn_genes);
    target_overlap{i,2}=intersect(targets{i},all_segregating_genes);
    
    temp_table=table([length(target_overlap{i,1});length(qtn_genes)-length(target_overlap{i,1})],...
        [length(target_overlap{i,2});length(all_segregating_genes)-length(target_overlap{i,2})],...
        'VariableNames',{'interacting','all other'},'RowNames',{'client','not client'});
    [h,p,stats]=fishertest(temp_table);
    
    qtn_pval(i)=p*length(kinases);
    
    
end

n_mat=cellfun(@length,target_overlap);

for i=1:length(kinases)
        
    f_mat(i,1)=n_mat(i,1)/length(qtn_genes);
    f_mat(i,2)=n_mat(i,2)/length(all_segregating_genes);  
     
end

qtn_odds_ratio=f_mat(:,1)./f_mat(:,2);

q_thresh=0.05;


to_use=qtn_pval<q_thresh;
q_to_use=qtn_pval(to_use);
kinase_to_use=kinases(to_use);


v1=qtn_odds_ratio(to_use);
v2=q_to_use;
v3=kinase_to_use;

[~,sort_idx]=sort(v1,'descend');

odds_for_merge=[];
q_for_merge=[];
labels_for_merge=[];


odds_for_merge=[odds_for_merge; v1(sort_idx)];
q_for_merge=[q_for_merge v2(sort_idx)];
labels_for_merge=[labels_for_merge; v3(sort_idx)];




%TF targets from huge ChIPseq experiment
%https://www.nature.com/articles/s41586-021-03314-8#MOESM4

[num,txt]=xlsread([dependency_directory '41586_2021_3314_MOESM4_ESM.xlsx'],6);

v_tf=txt(2,4:374);

v_tf_genes=txt(9:5872,2);

v_tf_type=txt(6,4:374);

clear txt

tf_mat=num(9:5872,4:374);

clear num

to_use=ismember(v_tf_type,'TF');

v_tf=v_tf(to_use);
tf_mat=tf_mat(:,to_use);


for i=1:length(v_tf)
    
    tf_targets{i}=v_tf_genes(logical(tf_mat(:,i)));
    
end

for i=1:length(v_tf)
    
    tf_target_overlap{i,1}=intersect(tf_targets{i},qtn_genes);
    tf_target_overlap{i,2}=intersect(tf_targets{i},all_segregating_genes);
    
    temp_table=table([length(tf_target_overlap{i,1});length(qtn_genes)-length(tf_target_overlap{i,1})],...
        [length(tf_target_overlap{i,2});length(all_segregating_genes)-length(tf_target_overlap{i,2})],...
        'VariableNames',{'interacting','all other'},'RowNames',{'client','not client'});
    [h,p,stats]=fishertest(temp_table);
    
    tf_qtn_pval(i)=p*length(v_tf);
    
    
end


n_mat_tf=cellfun(@length,tf_target_overlap);

for i=1:length(v_tf)
        
    f_mat_tf(i,1)=n_mat_tf(i,1)/length(qtn_genes);
    f_mat_tf(i,2)=n_mat_tf(i,2)/length(all_segregating_genes);  
     
end

tf_qtn_odds_ratio=f_mat_tf(:,1)./f_mat_tf(:,2);


q_thresh=0.05;
to_use=tf_qtn_pval<q_thresh;
q_to_use=tf_qtn_pval(to_use);
tf_to_use=v_tf(to_use);



%sort by descending OR for merged fig
v1=tf_qtn_odds_ratio(to_use);
v2=q_to_use;
v3=tf_to_use';

[~,sort_idx]=sort(v1,'descend');

odds_for_merge=[odds_for_merge; v1(sort_idx)];
q_for_merge=[q_for_merge v2(sort_idx)];
labels_for_merge=[labels_for_merge; v3(sort_idx)];



%also add enriched e3s
load([dependency_directory 'biogrid_data_physical.mat'])

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
    calculate_fraction_interactors(1,input_genes,dependency_directory,output_directory);

%Hsp90/70 machinery
%see https://journals.asm.org/doi/10.1128/mmbr.05018-11
list_to_use=4;
chaperone_input=readtable([dependency_directory 'chaperone_lists.csv']);
chaperone_query=table2array(chaperone_input(:,list_to_use));

chaperone_query(cellfun(@isempty,chaperone_query))=[];

for i=1:length(chaperone_query)
    chaperone_names{i}=all_labels{ismember(all_genes,chaperone_query{i})};
end

for i=1:length(chaperone_query)
    temp_idx=ismember(all_genes,chaperone_query{i});
    v_temp1(i)=relative_mat(1,temp_idx);
    v_temp2(i)=interactor_pval(temp_idx);
end

q_thresh=0.05;
to_use=find(v_temp2<q_thresh);
%sort by descending OR for merged fig
v1=v_temp1(to_use);
v2=v_temp2(to_use);
v3=chaperone_names(to_use);

[~,sort_idx]=sort(v1,'descend');



odds_for_merge=[odds_for_merge; v1(sort_idx)'];
q_for_merge=[q_for_merge v2(sort_idx)];
labels_for_merge=[labels_for_merge; v3(sort_idx)'];





hold on
bar(odds_for_merge)
%ylim([0 3])
set(gca,'YScale','log')
ylim([0.5 3])
plot(xlim,[1 1],':r')
ylabel('QTN odds ratio')
xticks(1:length(odds_for_merge))
xtickangle(45)
xticklabels(labels_for_merge)
for i=1:length(q_for_merge)
    if q_for_merge(i)<q_thresh
        h=text(i,2.2,num2str(q_for_merge(i)));
        set(h,'Rotation',45);
    end
end




end

