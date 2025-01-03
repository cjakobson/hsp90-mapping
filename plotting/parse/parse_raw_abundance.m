%parse strain names and raw abundance data
function [input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,...
    orf_names,strain_index]=parse_raw_abundance(dependency_directory,output_directory)

    input_data=readtable([dependency_directory '211031_SegregantProteomicsData_DetectionThreshold80_genes_ORF.tsv'],...
        'FileType','text');
    
    orf_names=input_data.Protein_Group;
    
    input_mat=table2array(input_data(:,2:end));
    %discard proteins with large number of missing values
    input_mat(sum(isnan(input_mat),2)>200,:)=nan;
    
    m=1;
    for i=2:length(input_data.Properties.VariableNames)

        temp_str=strsplit(input_data.Properties.VariableNames{i},'__');

        if length(temp_str)==3
            strain_names{m}=temp_str{1};
            ydj_names{m}=temp_str{2};
            strain_plate{m}=temp_str{3};
            strain_merge_idx{m}=erase(strain_plate{m},'_');
        else        
            temp_str2=strsplit(temp_str{1},'_');
            strain_names{m}=temp_str2{1};
            strain_index(m)=str2num(temp_str2{3});
        end

        m=m+1;

    end

    strain_merge_idx(cellfun(@isempty,strain_merge_idx))=[];

    rm_idx=find(ismember(strain_names,'RM11'));
    yjm_idx=find(ismember(strain_names,'YJM975'));
    f6_idx=find(ismember(strain_names,'F6'));
    
    sgrp_idx=find(~ismember(strain_names,'RM11').*...
        ~ismember(strain_names,'YJM975').*...
        ~ismember(strain_names,'F6'));

end

