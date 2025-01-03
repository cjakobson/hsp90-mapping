function [] = plot_interactors_abundance(input_data_switch,dependency_directory,output_directory)

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

if input_data_switch==0
    load([dependency_directory 'biogrid_data.mat'])
elseif input_data_switch==1
    load([dependency_directory 'biogrid_data_physical.mat'])
    interaction_mat=interaction_mat_physical;
elseif input_data_switch==2
    load([dependency_directory 'biogrid_data_genetic.mat'])
    interaction_mat=interaction_mat_genetic;
end



abundance_data=readtable([dependency_directory '4932/4932-WHOLE_ORGANISM-integrated.txt']);

n_interactors=nan(height(abundance_data),1);
for i=1:height(abundance_data)
    
    temp_str=strsplit(abundance_data.x_string_external_id{i},'.');
    
    abundance_orf{i}=temp_str{2};
    
    %look up number of interactors
    biogrid_idx=find(ismember(all_genes,abundance_orf{i}));
    
    if ~isempty(biogrid_idx)
        
        n_interactors(i)=sum(interaction_mat(biogrid_idx,:));
        
    end
    
end



v1=abundance_data.abundance;
v2=n_interactors;

scatter(v1,v2,10,'k','filled')
axis square
set(gca,'XScale','log')
set(gca,'YScale','log')
ylabel('# of interactors')
xlabel('est. protein abundance')
[r p]=corr(v1,v2,'rows','complete');
text(1e3,10e3,['r = ' num2str(r)])
text(1e3,7e3,['p = ' num2str(p)])




end


