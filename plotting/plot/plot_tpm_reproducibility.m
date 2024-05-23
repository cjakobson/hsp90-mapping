function [] = plot_tpm_reproducibility(dependency_directory,output_directory)


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



strains_to_plot={'RM11','YJM975'};

m=1;
for i=1:length(strains_to_plot)
    
    strain_idx=ismember(sample_info.strain,strains_to_plot{i});
    
    no_rad_idx=ismember(sample_info.condition,'untreated');
    rad_idx=ismember(sample_info.condition,'radicicol');
    
    idx1=logical(strain_idx.*no_rad_idx);
    temp_mat=tpm_mat(:,idx1);
    
    for j=1:3
        
        for k=(j+1):3
            
            v1=temp_mat(:,j);
            v2=temp_mat(:,k);
            
            subplot(2,6,m)
            hold on
            scatter(v1,v2,10,'k','filled','MarkerFaceAlpha',0.5)
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            xlim([1 1e6])
            ylim(xlim)
            axis square
            
            [r p]=corr(log10(v1),log10(v2),'rows','complete');
            text(10,1e5,num2str(r))
            
            m=m+1;
            
        end
        
    end
    
    
    idx2=logical(strain_idx.*rad_idx);
    temp_mat=tpm_mat(:,idx2);
    
    for j=1:3
        
        for k=(j+1):3
            
            v1=temp_mat(:,j);
            v2=temp_mat(:,k);
            
            subplot(2,6,m)
            hold on
            scatter(v1,v2,10,'k','filled','MarkerFaceAlpha',0.5)
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            xlim([1 1e6])
            ylim(xlim)
            axis square
            xlabel('TPM')
            ylabel('TPM')
            
            [r p]=corr(log10(v1),log10(v2),'rows','complete');
            text(10,1e5,num2str(r))
            
            m=m+1;
            
        end
        
    end
    

end



end

