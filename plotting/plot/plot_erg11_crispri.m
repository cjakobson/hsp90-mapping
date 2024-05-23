function [] = plot_erg11_crispri(plot_offset,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



exp_name='20230914sigmoid';

condition_names={'flc+aTc2','flc+aTc4','flc+aTc8','flc+aTc15','flc+aTc31',...
    'flc+aTc62','flc+aTc125','flc+aTc250','flc+aTc500','flc+aTc1000'...
    'flc+aTc2+rad','flc+aTc4+rad','flc+aTc8+rad','flc+aTc15+rad','flc+aTc31+rad',...
    'flc+aTc62+rad','flc+aTc125+rad','flc+aTc250+rad','flc+aTc500+rad','flc+aTc1000+rad'...
    'teb+aTc2','teb+aTc4','teb+aTc8','teb+aTc15','teb+aTc31',...
    'teb+aTc62','teb+aTc125','teb+aTc250','teb+aTc500','teb+aTc1000'...
    'teb+aTc2+rad','teb+aTc4+rad','teb+aTc8+rad','teb+aTc15+rad','teb+aTc31+rad',...
    'teb+aTc62+rad','teb+aTc125+rad','teb+aTc250+rad','teb+aTc500+rad','teb+aTc1000+rad'};
    

strain_names={'YJM975 no guide','RM11 no guide','YJM975 antiERG11','RM11 antiERG11'};

v_time={'48h'};


n_images=14;
n_plates=40;

plate_names={'a','b','c'};


m=1;
%get sga data
mat_to_process=nan(1536,n_plates);
for i=1:n_images
    
    for j=1:length(plate_names)
        
        if length(num2str(i))==1
            temp_name=[exp_name '00' num2str(i) plate_names{j} '.jpg.dat'];
        elseif length(num2str(i))==2
            temp_name=[exp_name '0' num2str(i) plate_names{j} '.jpg.dat'];
        end    

        if exist([dependency_directory 'sigmoid/gitter/' temp_name])

            sga_mat{m}=readtable([dependency_directory 'sigmoid/gitter/' temp_name]);

            mat_to_process(:,m)=sga_mat{m}.size;

        end
        
        m=m+1;
        
    end
    
end


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

reorder_mat((1:384),:)=mat_to_process(a1idx,:);
reorder_mat(((384+1):(2*384)),:)=mat_to_process(a2idx,:);
reorder_mat(((2*384+1):(3*384)),:)=mat_to_process(b1idx,:);
reorder_mat(((3*384+1):(4*384)),:)=mat_to_process(b2idx,:);
        
reorder_mat(reorder_mat<50)=nan;




v_mean_all=[];
v_sem_all=[];
to_plot_all=[];

for i=1:length(v_time)
    
    m=1;
    
    for j=1:length(condition_names)
        
        columnOffset=(i-1)*length(condition_names);
        idxToUse=columnOffset+j;
        
        clear to_plot
        
        for k=1:length(strain_names)
            
            to_plot{k}=reorder_mat(384*(k-1)+(1:384),idxToUse);
            
        end
        
        %store raw values before normalization
        to_plot_all=[to_plot_all to_plot];
        
        %normalize to no edit
        v_temp=to_plot{1};
        for ii=[1 3]%1:length(to_plot)
            to_plot{ii}=to_plot{ii}./v_temp;
        end
        v_temp=to_plot{2};
        for ii=[2 4]%1:length(to_plot)
            to_plot{ii}=to_plot{ii}./v_temp;
        end
        
        
        for ii=1:length(to_plot)
            
            to_plot{ii}(isinf(to_plot{ii}))=nan;
            
            v_mean(ii)=mean(to_plot{ii},'omitnan');
            v_sem(ii)=std(to_plot{ii},[],'omitnan')./sqrt(length(to_plot{ii}));
        end
        
        
        v_mean_all=[v_mean_all v_mean];
        v_sem_all=[v_sem_all v_sem];
        
                
    end
    
    
end



%separate out RM and YJM
temp_labels1={'YJM975 antiERG11','RM11 antiERG11'};
temp_labels2={'fluconazole','tebuconazole'};
temp_labels3={'2ng/uL','4ng/uL','8ng/uL','15ng/uL','31ng/uL',...
    '62ng/uL','125ng/uL','250ng/uL','500ng/uL','1000ng/uL'};

m=1;

for i=1:length(temp_labels1)
    
    v_mean_to_plot=v_mean_all((2+i):4:length(v_mean_all));
    v_sem_to_plot=v_sem_all((2+i):4:length(v_sem_all));
    v_to_plot_all=to_plot_all((i):2:length(v_sem_all));
    
    for j=1:length(temp_labels2)
        
        temp_idx=20*(j-1)+(1:20);
        
        v1=v_mean_to_plot(temp_idx);
        v2=v_sem_to_plot(temp_idx);
        
        temp_idx2=20*(j-1)+(1:40);
        v3=v_to_plot_all(temp_idx2);
        
        
        if (i==1)&&(j==1)
            %also plot boxplot of raw values
            v3=[v3(11:20) v3(31:40)];

            %interleave -/+rad
            for ii=1:5
                v4(4*(5-ii)+(1:2))=v3(2*(ii-1)+(1:2));
                v4(4*(5-ii)+(3:4))=v3(2*(ii-1)+(1:2)+10);
            end

        end
        
        m=m+1;
            
    end
    
end




m=1;

%subplots of v4 for figures
for i=1:4
    
    to_sub_plot{i}=v4(i:4:end);
    
    for j=1:length(to_sub_plot{i})
        
        if i<=2
            v_temp=to_sub_plot{1}{j};
        else
            v_temp=to_sub_plot{3}{j};
        end
        
        v_sub_mean{i}(j)=mean(to_sub_plot{i}{j}./v_temp,'omitnan');
        v_sub_std{i}(j)=std(to_sub_plot{i}{j}./v_temp,[],'omitnan')./sqrt(length(to_sub_plot{i}{j}));
        
    end
    
end

%no rad then rad
n=1;
for i=[1 3]
    
    m=1;
    for j=1:length(to_sub_plot{i})

        to_sub_plot2{i}(m)=to_sub_plot{i}(j);
        m=m+1;
        
        to_sub_plot2{i}(m)=to_sub_plot{i+1}(j);
        m=m+1;
        
    end
    
    subplot(2,4,n+plot_offset)
    hold on
    easy_box(to_sub_plot2{i})
    ylim([0 350])
    plot(xlim,[median(to_sub_plot2{i}{end}) median(to_sub_plot2{i}{end})],':r')
    ylabel('spot size')
    n=n+1;

    
end




end


