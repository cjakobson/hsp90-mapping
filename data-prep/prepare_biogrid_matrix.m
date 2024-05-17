function [] = prepare_biogrid_matrix(dependency_directory,output_directory)

tic

biogrid_data=readtable([dependency_directory 'BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-4.4.207.tab3.txt']);

all_genes=unique([biogrid_data.SystematicNameInteractorA;...
    biogrid_data.SystematicNameInteractorB]);
%get rid of non-systematic things
all_genes=all_genes(462:6350);


all_labels=cell(size(all_genes));
for i=1:length(all_genes)

    if mod(i,100)==0
        i
    end
            
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
    
    if mod(i,100)==0
        i
    end
    
    query_idx=logical(ismember(biogrid_data.SystematicNameInteractorA,all_genes{i})+...
        ismember(biogrid_data.SystematicNameInteractorB,all_genes{i}));
    
    query_interactors=unique([biogrid_data.SystematicNameInteractorA(query_idx);...
        biogrid_data.SystematicNameInteractorB(query_idx)]);
        
    temp_idx=ismember(all_genes,query_interactors);
    
    interaction_mat(i,temp_idx)=1;
    
end

%imagesc(interaction_mat)

%output gene lists and matrix
save([output_directory 'biogrid_data.mat'],'all_genes','all_labels','interaction_mat')

toc


end

