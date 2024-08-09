function [] = plot_1K_effect(locus_to_plot,condition_to_plot,time_to_plot,dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;





exp_name={'20221030phenotyping','20221031phenotyping','20221101phenotyping','20221102phenotyping'};
exp_time={'24h','48h','72h','96h'};

n_plates=102;


input_data=[];
for i=1:length(exp_name)
    temp_table=readtable([dependency_directory exp_name{i} 'data.csv']);
    input_data=[input_data table2array(temp_table)];
end

[num txt]=xlsread([dependency_directory '20211209 radicicol rescreen plate key_FINAL.xlsx']);



temp_conditions1=txt(2:13,2);


full_mat=[];
for i=1:(n_plates*length(exp_time))
        
    
    temp_mat{i}=input_data(:,i);
    
    %rearrange
    %reorganize to 384
    a1idxBase=1:2:48;
    a2idxBase=2:2:48;
    b1idxBase=49:2:96;
    b2idxBase=50:2:96;

    a1idx=[];
    a2idx=[];
    b1idx=[];
    b2idx=[];

    for k=1:16

        a1idx=[a1idx 96*(k-1)+a1idxBase];
        a2idx=[a2idx 96*(k-1)+a2idxBase];
        b1idx=[b1idx 96*(k-1)+b1idxBase];
        b2idx=[b2idx 96*(k-1)+b2idxBase];

    end

    reorder_mat{i}=nan(size(temp_mat{i}));

    reorder_mat{i}((1:384))=temp_mat{i}(a1idx);
    reorder_mat{i}(((384+1):(2*384)))=temp_mat{i}(a2idx);
    reorder_mat{i}(((2*384+1):(3*384)))=temp_mat{i}(b1idx);
    reorder_mat{i}(((3*384+1):(4*384)))=temp_mat{i}(b2idx);

    full_mat=[full_mat reorder_mat{i}];
    
end

full_mat(full_mat<25)=nan;
full_mat(full_mat>2000)=nan;


%trim to just 1K
plate_offset=102;
plate_range=1:48;

to_keep=[];
for i=1:length(exp_time)
    to_keep=[to_keep (i-1)*plate_offset+plate_range];
end

full_mat=full_mat(:,to_keep);

n_plates=48;

%get strain info
strain_info=readtable([dependency_directory '1011 Genomes_Sace_strains_matrix_positions_384.xlsx']);
m=1;
for i=1:3   %plates
    
    temp_idx1=ismember(strain_info.Matrix_384,['M' num2str(i)]);
    
    for j=1:16  %rows
        
        temp_idx2=ismember(strain_info.row_384,j);
        
        for k=1:24  %columns
            
            temp_idx3=ismember(strain_info.col_384,k);
            
            temp_idx=find(temp_idx1.*temp_idx2.*temp_idx3);
            
            if length(temp_idx)>0
                
                name1{m}=strain_info.Strain_Name_4SRA{temp_idx};
                name2{m}=strain_info.Standardized_name{temp_idx};
                
            else
                
                name1{m}='NA';
                name2{m}='NA';
                
            end
            
            m=m+1;
            
        end
        
    end
        
end


%should recast this to be deltaGrowth -/+ rad
%calculate change in growth for each replicate
for i=1:(n_plates*length(exp_time)/4)
    
    idx1=(i-1)*4+1;
    idx2=(i-1)*4+3;
    
    no_rad_mean_mat(:,i)=mean([full_mat(:,idx1) full_mat(:,idx2)],2);
    
    idx1=(i-1)*4+2;
    idx2=(i-1)*4+4;
    
    rad_mean_mat(:,i)=mean([full_mat(:,idx1) full_mat(:,idx2)],2);
    
end


   

condition_to_use=find(ismember(temp_conditions1,condition_to_plot));
time_to_use=find(ismember(exp_time,time_to_plot));

%minor allele frequency from 1k genomes
load([dependency_directory '1002data.mat'])

temp_het_mat=minGenotype~=minGenotype2;

minGenotype(temp_het_mat)=-1;

for i=1:length(locus_to_plot)

    altIdx=logical((minGenotype(locus_to_plot(i),:)==1)+(minGenotype(locus_to_plot(i),:)==-1));
    
    altStrains{i}=strainString(altIdx);

end


%normalize to plot
for j=time_to_use
    
    column_offset=n_plates/4*(j-1);
    
    for i=condition_to_use

        condition_idx=i;
        hold on
        to_plot{2}=no_rad_mean_mat(1:1152,condition_idx+column_offset);

        altIdx=find(ismember(name2,altStrains{1}));

        to_plot{1}=to_plot{2}(altIdx);
        to_plot{2}(altIdx)=[];
        
        to_plot{4}=rad_mean_mat(1:1152,condition_idx+column_offset);
        
        to_plot{3}=to_plot{4}(altIdx);
        to_plot{4}(altIdx)=[];
        
        clear v_mean v_sem
        for k=1:length(to_plot)
            v_mean(k)=median(to_plot{k},'omitnan');
            v_sem(k)=std(to_plot{k},[],'omitnan')./sqrt(length(to_plot{k}));
        end
        
        v_sem(1:2)=v_sem(1:2)./v_mean(1);
        v_mean(1:2)=v_mean(1:2)./v_mean(1);
        
        v_sem(3:4)=v_sem(3:4)./v_mean(3);
        v_mean(3:4)=v_mean(3:4)./v_mean(3);
            
        
        bar(v_mean)
        errorbar(1:length(v_mean),v_mean,v_sem,'k.')

        ylabel('relative growth')
        title([condition_to_plot ' ' time_to_plot])
        temp_labels={'ref-rad','alt-rad','ref+rad','alt+rad'};
        xticks(1:length(temp_labels))
        xticklabels(temp_labels)
        xtickangle(45)
        [h p]=ttest2(to_plot{2},to_plot{4});
        text(3,1.7,num2str(p))
        xlim([0 5])

        for k=1:length(to_plot)
            text(k,0.1,num2str(length(to_plot{k})))
        end
        

    end
    

end




end



