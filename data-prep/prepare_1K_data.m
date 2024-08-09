function [] = prepare_1K_data(dependency_directory,output_directory)

tic


load([dependency_directory '1002_data_all_variants.mat'])

gene(cellfun(@isempty,gene))={'NA'};
type(cellfun(@isempty,type))={'NA'};


all_genes=unique(gene);
all_genes=all_genes(2:end);




%get number of mutations and allele frequency spectrum
mis_n_muts=nan(length(all_genes),1);
syn_n_muts=nan(length(all_genes),1);
up_n_muts=nan(length(all_genes),1);
down_n_muts=nan(length(all_genes),1);


mis_af_gene=cell(length(all_genes),1);
syn_af_gene=cell(length(all_genes),1);
up_af_gene=cell(length(all_genes),1);
down_af_gene=cell(length(all_genes),1);


mis_mean_af_gene=nan(length(all_genes),1);
syn_mean_af_gene=nan(length(all_genes),1);
up_mean_af_gene=nan(length(all_genes),1);
down_mean_af_gene=nan(length(all_genes),1);


mis_idx=ismember(type,'missense_variant');
syn_idx=ismember(type,'synonymous_variant');
reg_idx=logical(ismember(type,'downstream_gene_variant')+...
    ismember(type,'upstream_gene_variant'));

%normalize to gene length
gene_info=readtable([dependency_directory 'TableS4.xls']);

%dist_thresh=500;

for i=1:length(all_genes)

    if mod(i,100)==0

        i

    end

    temp_gene_idx=ismember(gene,all_genes{i});

    temp_mis_idx=ismember(type(temp_gene_idx),'missense_variant');
    temp_syn_idx=ismember(type(temp_gene_idx),'synonymous_variant');

    temp_gene_info_idx=find(ismember(gene_info.Name,all_genes{i}));

    if ~isempty(temp_gene_info_idx)
        
        temp_length=abs(gene_info.SGD_Start(temp_gene_info_idx)-...
            gene_info.SGD_End(temp_gene_info_idx));
    
        mis_n_muts(i)=sum(temp_mis_idx)/temp_length;

        %also incorporate allele frequency
        mis_af_gene{i}=af(logical(temp_gene_idx.*mis_idx));
        mis_mean_af_gene(i)=mean(mis_af_gene{i},'omitnan');
        
        
        %also synonymous
        syn_n_muts(i)=sum(temp_syn_idx)/temp_length;
        
        syn_af_gene{i}=af(logical(temp_gene_idx.*syn_idx));
        syn_mean_af_gene(i)=mean(syn_af_gene{i},'omitnan');
        
        %calculate distance
        
        %just consider upstream region
        temp_chr_idx=ismember(chr,gene_info.Chrom{temp_gene_info_idx}(4:end));
        
        temp_start_pos=gene_info.SGD_Start(temp_gene_info_idx);
        temp_end_pos=gene_info.SGD_End(temp_gene_info_idx);
        
        if temp_start_pos<temp_end_pos  %watson
            
            temp_up_dist=temp_start_pos-pos;
            temp_up_dist(temp_up_dist<0)=nan;
            
            temp_down_dist=pos-temp_end_pos;
            temp_down_dist(temp_down_dist<0)=nan;
            
            
            temp_up_gene_idx=max([temp_gene_info_idx-1,1]);
            
            temp_down_gene_idx=min([temp_gene_info_idx+1,height(gene_info)]);
            

            temp_up_next_gene_pos=max([gene_info.SGD_Start(temp_up_gene_idx),...
                gene_info.SGD_End(temp_up_gene_idx)]);
            temp_up_next_gene_dist=abs(temp_up_next_gene_pos-temp_start_pos);
            
            temp_down_next_gene_pos=min([gene_info.SGD_Start(temp_down_gene_idx),...
                gene_info.SGD_End(temp_down_gene_idx)]);
            temp_down_next_gene_dist=abs(temp_down_next_gene_pos-temp_end_pos);
            
            
            temp_dist_idx=temp_up_dist<temp_up_next_gene_dist;
            
            up_n_muts(i)=sum(reg_idx.*temp_dist_idx.*temp_chr_idx)/...
                temp_up_next_gene_dist;
            
            up_af_gene{i}=af(logical(reg_idx.*temp_dist_idx.*temp_chr_idx));
            up_mean_af_gene(i)=mean(up_af_gene{i},'omitnan');
            
            
            temp_dist_idx=temp_down_dist<temp_down_next_gene_dist;
            
            down_n_muts(i)=sum(reg_idx.*temp_dist_idx.*temp_chr_idx)/...
                temp_down_next_gene_dist;
            
            down_af_gene{i}=af(logical(reg_idx.*temp_dist_idx.*temp_chr_idx));
            down_mean_af_gene(i)=mean(down_af_gene{i},'omitnan');
            
            
        else    %crick
            
            temp_up_dist=temp_end_pos-pos;
            temp_up_dist(temp_up_dist<0)=nan;
            
            temp_down_dist=pos-temp_start_pos;
            temp_down_dist(temp_down_dist<0)=nan;
            
            
            temp_up_gene_idx=max([temp_gene_info_idx-1,1]);
            
            temp_down_gene_idx=min([temp_gene_info_idx+1,height(gene_info)]);
            
        
            temp_up_next_gene_pos=max([gene_info.SGD_Start(temp_up_gene_idx),...
                gene_info.SGD_End(temp_up_gene_idx)]);
            temp_up_next_gene_dist=abs(temp_up_next_gene_pos-temp_end_pos);
            
            temp_down_next_gene_pos=min([gene_info.SGD_Start(temp_down_gene_idx),...
                gene_info.SGD_End(temp_down_gene_idx)]);
            temp_down_next_gene_dist=abs(temp_down_next_gene_pos-temp_start_pos);
            
            
            temp_dist_idx=temp_down_dist<temp_down_next_gene_dist;
            
            up_n_muts(i)=sum(reg_idx.*temp_dist_idx.*temp_chr_idx)/...
                temp_down_next_gene_dist;
            
            up_af_gene{i}=af(logical(reg_idx.*temp_dist_idx.*temp_chr_idx));
            up_mean_af_gene(i)=mean(up_af_gene{i},'omitnan');
            
            
            temp_dist_idx=temp_up_dist<temp_up_next_gene_dist;
            
            down_n_muts(i)=sum(reg_idx.*temp_dist_idx.*temp_chr_idx)/...
                temp_up_next_gene_dist;
            
            down_af_gene{i}=af(logical(reg_idx.*temp_dist_idx.*temp_chr_idx));
            down_mean_af_gene(i)=mean(down_af_gene{i},'omitnan');
           
        end
        
        
        
    end

end


%exclude those with no mutations
mis_n_muts(mis_n_muts==0)=nan;
syn_n_muts(syn_n_muts==0)=nan;
up_n_muts(up_n_muts==0)=nan;
down_n_muts(down_n_muts==0)=nan;

%take mean of up and down as regulatory

reg_n_muts=mean([up_n_muts down_n_muts],2);

%save as table

to_output=table(all_genes',mis_n_muts,reg_n_muts);
writetable(to_output,[output_directory '1K_frequencies.csv'])

toc


end

