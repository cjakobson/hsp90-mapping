function [] = prepare_biogrid_matrix(dependency_directory,output_directory)

tic

biogrid_data=readtable([dependency_directory 'BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.207.tab3.txt']);

all_genes=unique([biogrid_data.SystematicNameInteractorA;...
    biogrid_data.SystematicNameInteractorB]);
%get rid of non-systematic things
all_genes=all_genes(462:6350);


all_labels=cell(size(all_genes));
for i=1:length(all_genes)
            
    v_temp=ismember(biogrid_data.SystematicNameInteractorA,all_genes{i});
    if sum(v_temp)>0
        
        v_temp=biogrid_data.OfficialSymbolInteractorA(v_temp);
        all_labels{i}=v_temp{1};
        
    else
        
        all_labels{i}=' ';
        
    end
        
end


%build matrix of interactions
interaction_mat=zeros(length(all_genes));

for i=1:length(all_genes)
    
    query_idx=logical(ismember(biogrid_data.SystematicNameInteractorA,all_genes{i})+...
        ismember(biogrid_data.SystematicNameInteractorB,all_genes{i}));
    
    query_interactors=unique([biogrid_data.SystematicNameInteractorA(query_idx);...
        biogrid_data.SystematicNameInteractorB(query_idx)]);
        
    temp_idx=ismember(all_genes,query_interactors);
    
    interaction_mat(i,temp_idx)=1;
    
end

%output gene lists and matrix
save([output_directory 'biogrid_data.mat'],'all_genes','all_labels','interaction_mat')


%repeat for physical and genetic interactions separately
biogrid_data_physical=biogrid_data(ismember(biogrid_data.ExperimentalSystemType,'physical'),:);
biogrid_data_genetic=biogrid_data(ismember(biogrid_data.ExperimentalSystemType,'genetic'),:);

interaction_mat_physical=zeros(length(all_genes));
interaction_mat_genetic=zeros(length(all_genes));

for i=1:length(all_genes)
    
    query_idx_physical=logical(ismember(biogrid_data_physical.SystematicNameInteractorA,all_genes{i})+...
        ismember(biogrid_data_physical.SystematicNameInteractorB,all_genes{i}));
    
    query_interactors_physical=unique([biogrid_data_physical.SystematicNameInteractorA(query_idx_physical);...
        biogrid_data_physical.SystematicNameInteractorB(query_idx_physical)]);
        
    temp_idx=ismember(all_genes,query_interactors_physical);

    if sum(temp_idx)>0
    
        interaction_mat_physical(i,temp_idx)=1;
        
    end

    
    query_idx_genetic=logical(ismember(biogrid_data_genetic.SystematicNameInteractorA,all_genes{i})+...
        ismember(biogrid_data_genetic.SystematicNameInteractorB,all_genes{i}));
    
    query_interactors_genetic=unique([biogrid_data_genetic.SystematicNameInteractorA(query_idx_genetic);...
        biogrid_data_genetic.SystematicNameInteractorB(query_idx_genetic)]);
        
    temp_idx=ismember(all_genes,query_interactors_genetic);

    if sum(temp_idx)>0
    
        interaction_mat_genetic(i,temp_idx)=1;
        
    end
end


save([output_directory 'biogrid_data_physical.mat'],'all_genes','all_labels',...
    'interaction_mat_physical')
save([output_directory 'biogrid_data_genetic.mat'],'all_genes','all_labels',...
    'interaction_mat_genetic')



toc


end

