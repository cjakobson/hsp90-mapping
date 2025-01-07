function [] = plot_heritability_explained_norm(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

mapping_input=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);


exp_name={'20221030phenotyping','20221031phenotyping','20221101phenotyping','20221102phenotyping'};
exp_time={'24h','48h','72h','96h'};

n_plates=102;

input_data=[];
for i=1:length(exp_name)
    tempTable=readtable([dependency_directory exp_name{i} 'data.csv']);
    input_data=[input_data table2array(tempTable)];
end

[num txt]=xlsread([dependency_directory '20211209 radicicol rescreen plate key_FINAL.xlsx']);



temp_conditions1=txt(2:13,2);
temp_conditions2=num(1:12,1);
temp_conditions3=txt(2:13,4);

for i=1:length(temp_conditions1)
    temp_conditions{i}=[temp_conditions1{i} '_' num2str(temp_conditions2(i))...
        temp_conditions3{i}];
end


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


%trim to just heritability/DDO part for -/+rad
plateOffset=102;
plateRange=49:96;

toKeep=[];
for i=1:length(exp_time)
    toKeep=[toKeep (i-1)*plateOffset+plateRange];
end

full_mat=full_mat(:,toKeep);

[~, n_plates]=size(full_mat);
for i=1:n_plates
    
    v_temp=full_mat(:,i);
    
    z_mat(:,i)=(v_temp-mean(v_temp,'omitnan'))./std(v_temp,[],'omitnan');
    
end


%order is-rad rep 1 +rad rep1 -rad rep 2 +rad rep 2 for 12 conditions
for i=1:2*length(exp_time)*length(temp_conditions1)
    
    v1=z_mat(:,2*(i-1)+1);
    v2=z_mat(:,2*(i-1)+2);
    
    delta_mat(:,i)=v2-v1;
    
end



m=1;
for i=1:length(exp_time)
    
    for j=1:length(temp_conditions1)
        
        temp_labels{m}=[exp_time{i} ' ' temp_conditions1{j}];
        m=m+1;
        
    end
    
end

for i=1:length(exp_time)*length(temp_conditions1)
    
    delta1=delta_mat(:,2*(i-1)+1);
    delta2=delta_mat(:,2*(i-1)+2);
    
    if (sum(~isnan(delta1))>0)&&(sum(~isnan(delta2))>0)
        temp_table=table(delta1,delta2);
        temp_lme=fitlme(temp_table,'delta1~delta2');
    end    
    h_squared(i)=1-temp_lme.MSE;
        
end

geo_idx=mapping_input.index==0;

m=1;
for i=1:length(exp_time)
    
    time_idx=ismember(mapping_input.time,exp_time{i});
    
    for j=1:length(temp_conditions)
        
        condition_idx=ismember(mapping_input.condition,temp_conditions{j});
        
        tempIdx=logical(time_idx.*condition_idx);
        
        v_geo(m)=sum(mapping_input.varExp(logical(tempIdx.*geo_idx)));
        v_linear(m)=sum(mapping_input.varExp(logical(tempIdx.*~geo_idx)));
        
        v_remainder(m)=h_squared(m)-(v_geo(m)+v_linear(m));
        
        f_explained(m)=(v_geo(m)+v_linear(m))/h_squared(m);
        
        m=m+1;
        
    end
    
end

%sort(f_explained,'ascend')

%SDS didn't map well
to_remove=11:12:length(temp_labels);
to_remove=ismember(1:length(temp_labels),to_remove);

temp_mat=[v_geo(~to_remove); v_linear(~to_remove); v_remainder(~to_remove)]';


tempIdx=[1:10 12];

%plot all 4 time points?

m=1;
for i=4%1:length(exp_time)

    %tempMat2=temp_mat((end-10):end,:);
    tempMat2=temp_mat((11*(i-1)+1):(11*i),:);
    [~,sortIdx]=sort(tempMat2(:,2),'ascend');
    
    
    %norm to H2
    temp_mat3=tempMat2(sortIdx,1:3);
    v_temp=sum(temp_mat3,2);
    

    v1=temp_mat3(:,2)./v_temp;
    v2=sum(temp_mat3(:,1:2),2)./v_temp;
    %subplot(2,3,plot_offset+m)
    m=m+1;
    hold on
    %bar(tempMat2(sortIdx,:),'stacked')
    bar(v1,'stacked')
    max(v1)
    %max(v2)
    ylim([0 1])
    xticks(1:length(tempIdx))
    xtickangle(45)
    temp_labels=temp_conditions(tempIdx);
    xticklabels(temp_labels(sortIdx))
    title(exp_time{i})
    %legend({'Hsp90-dependent'})
    axis square
    ylabel('fraction of H^2 explained')

end


end


