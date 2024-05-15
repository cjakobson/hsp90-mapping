function [] = plot_heat_shock_rnas(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

sample_info=readtable([dependency_directory 'sampleInfo.txt']);

%ignore het
sample_info(7:12,:)=[];

filebase=[dependency_directory 'tpms/'];

%load data
for i=1:height(sample_info)
    
    temp_data=tdfread([filebase sample_info.sample{i} '/abundance.tsv']);
    tpm_mat(:,i)=temp_data.tpm;

end



%remove low-abundance transcripts
tpm_mat(tpm_mat<1)=nan;

genes=temp_data.target_id;

for i=1:length(genes)
    
    temp_str=strsplit(genes(i,:),'_');
    tpm_genes{i}=temp_str{1};
    
end



strains=unique(sample_info.strain);
conditions=unique(sample_info.condition);



%analyze spike-ins then exclude before PCA
erccNames=tpm_genes(6601:end);
erccData=tpm_mat(6601:end,:);

tpm_genes=tpm_genes(1:6600);
tpm_mat=tpm_mat(1:6600,:);


%plot some relevant chaperone genes explicitly (TPMs)
%see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4938784/
common_to_plot={'HSC82','HSP82','SSA1','SSA2','HSP78','HSP104','SIS1','BTN2'};
systematic_to_plot={'YMR186W','YPL240C','YAL005C','YLL024C',...
    'YDR258C','YLL026W','YNL007C','YGR142W'};

strains_to_plot={'YJM975','RM11'};

clear v1
clear v2
for i=1:length(strains_to_plot)
    
    strain_idx=ismember(sample_info.strain,strains_to_plot{i});
    
    no_rad_idx=ismember(sample_info.condition,'untreated');
    rad_idx=ismember(sample_info.condition,'radicicol');
    
    m=1;
    for j=1:length(systematic_to_plot)
        
        gene_idx=ismember(tpm_genes,systematic_to_plot{j});
        
        idx1=logical(strain_idx.*no_rad_idx);
        
        v1{m}=tpm_mat(gene_idx,idx1);
        v2(m)=mean(v1{m},'omitnan');
        temp_labels{m}=[common_to_plot{j} ' noRad'];
        m=m+1;
        
        idx2=logical(strain_idx.*rad_idx);
        
        v1{m}=tpm_mat(gene_idx,idx2);
        v2(m)=mean(v1{m},'omitnan');
        temp_labels{m}=[common_to_plot{j} ' rad'];
        m=m+1;
        
    end
    
    
    %FC
    v3=v2(1:2:end);
    v4=v2(2:2:end);
    
    v5=v1(1:2:end);
    v6=v1(2:2:end);
    
    subplot(2,4,i+2)
    hold on
    bar(v4./v3)
    for j=1:length(v5)
        scatter(ones(length(v5{j}),1)*j,v6{j}./v5{j},10,'k','filled')
    end
    xticks(1:length(common_to_plot))
    xticklabels(common_to_plot)
    xtickangle(45)
    ylabel('FC +rad/-rad')
    title(strains_to_plot{i})
    ylim([0.1 10])
    set(gca,'YScale','log')
    plot(xlim,[1 1],':k')

end



end

