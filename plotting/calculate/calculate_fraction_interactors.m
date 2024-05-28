function [overlap_mat,interactor_pval,relative_mat] = ...
    calculate_fraction_interactors(input_data_switch,input_genes,dependency_directory,output_directory)



if input_data_switch==0
    load([dependency_directory 'biogrid_data.mat'])
elseif input_data_switch==1
    load([dependency_directory 'biogrid_data_physical.mat'])
    interaction_mat=interaction_mat_physical;
elseif input_data_switch==2
    load([dependency_directory 'biogrid_data_genetic.mat'])
    interaction_mat=interaction_mat_genetic;
end


interaction_thresh=25;

overlap_mat=zeros(2,length(all_genes));
n_overlap=overlap_mat;
interactor_pval=nan(length(all_genes),1);
for j=1:length(all_genes)
    
    v_temp=logical(interaction_mat(j,:));
    
    if sum(v_temp)>interaction_thresh
        
        query_interactors=all_genes(v_temp);
        
        %count overlap if either associated QTN gene has overlap
        overlap_mat(1,j)=sum(logical(ismember(input_genes{1},query_interactors)+...
            ismember(input_genes{2},query_interactors)))./length(input_genes{1});
        n_overlap(1,j)=sum(logical(ismember(input_genes{1},query_interactors)+...
            ismember(input_genes{2},query_interactors)));
            
        %all other
        overlap_mat(2,j)=sum(ismember(input_genes{3},query_interactors))./length(input_genes{3});
        n_overlap(2,j)=sum(ismember(input_genes{3},query_interactors));
        
        %calculate p value for hsp90 vs all other segregating
        temp_table=table([n_overlap(1,j);length(input_genes{1})-n_overlap(1,j)],...
            [n_overlap(2,j);length(input_genes{3})-n_overlap(2,j)],...
            'VariableNames',{'interacting','all other'},'RowNames',{'client','not client'});
        [~,p,~]=fishertest(temp_table);
        interactor_pval(j)=p;
        
    end
        
end

%Bonferroni
interactor_pval=interactor_pval.*sum(~isnan(interactor_pval));

for i=1:2
    relative_mat(i,:)=overlap_mat(i,:)./overlap_mat(end,:);
end

end


