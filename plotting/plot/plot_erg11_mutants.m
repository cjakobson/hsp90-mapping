function [] = plot_erg11_mutants(condition_to_plot,dependency_directory,output_directory)




set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



%read scanner sgaData files


start_date=20230727;
start_time=12*60+24+13/60; 

plate_names={'a','b','c','d'};%,'e'};

n_scanners=3;

condition_names={'flc-rad','flc+rad','teb-rad','teb+rad'};

strain_names{1}={'6635_YJM975','6649_RM11','8281_YJM975_ERG11-1','8280_RM11_ERG11-1'};
strain_names{2}={'6635_YJM975','8281_YJM975_ERG11-1','8436_YJM975_ERG11-2','8437_YJM975_ERG11-1+2'};
strain_names{3}={'6649_RM11','8280_RM11_ERG11-1','8438_RM11_ERG11-2','8439_RM11_ERG11-1+2'};


%get filenames
to_read=dir([dependency_directory 'erg11rad/gitter']);


%get sga data
m=1;

for i=1:length(to_read)
    
    temp_name=to_read(i).name;

    if length(temp_name)>3

        if strcmp(temp_name((end-2):end),'dat')

            to_keep{m}=temp_name;

            %grab dates and times
            temp_str=strsplit(temp_name,'_');
            dates(m)=str2num(temp_str{1});
            times{m}=temp_str{end-2};

            %convert time to minutes and normalize dates
            temp_date=dates(m)-start_date;

            temp_hours=times{m}(1:2);
            temp_mins=times{m}(3:4);
            temp_secs=times{m}(5:6);


            time_mins(m)=str2num(temp_hours)*60+str2num(temp_mins)+...
                str2num(temp_secs)/60+24*60*temp_date-start_time;

            temp_str2=strsplit(temp_str{end},'.');

            scanner(m)=str2num(temp_str2{1}(1));
            plates{m}=temp_str2{1}(2:end);

            %get sga data
            sga_mat{m}=readtable([dependency_directory 'erg11rad/gitter/' temp_name]);

            mat_to_process(:,m)=sga_mat{m}.size;


            m=m+1;

        end

    end

end


m=1;
for k=1:n_scanners
    

    for i=1:length(plate_names)

        scanner_idx=scanner==k;
        plate_idx=ismember(plates,plate_names{i});
        time_idx=time_mins>0;

        [n_strains,~]=size(mat_to_process);

        idx_to_use=logical(scanner_idx.*plate_idx.*time_idx);

        temp_mat=mat_to_process(:,idx_to_use);

        %remove discontinuities
        temp_delta_mat=temp_mat(:,2:end)-temp_mat(:,1:(end-1));
        to_remove=temp_delta_mat<-50;
        to_remove=[zeros(1536,1) to_remove];
        to_remove=logical(to_remove);

        temp_mat(to_remove)=nan;
        
        growth_mat{m}=temp_mat;

        m=m+1;

    end

end



%rearrange everything to 384 first
for i=1:length(growth_mat)
    
    temp_mat=growth_mat{i};
    
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

    reorder_mat{i}=nan(size(temp_mat));

    reorder_mat{i}((1:384),:)=temp_mat(a1idx,:);
    reorder_mat{i}((384+1):(2*384),:)=temp_mat(a2idx,:);
    reorder_mat{i}((2*384+1):(3*384),:)=temp_mat(b1idx,:);
    reorder_mat{i}((3*384+1):(4*384),:)=temp_mat(b2idx,:);
    
end




fwhmIdx=48;





%p values for figure
v_mean_all=[];
v_sem_all=[];
to_plot_all=[];

m=1;
for i=1:length(strain_names)
    
    for j=1:length(condition_names)
        
        temp_mat=reorder_mat{m};
        
        v_to_use=temp_mat(:,fwhmIdx);
        
        clear to_plot
        for k=1:length(strain_names{i})
            
            temp_idx=384*(k-1)+(1:384);
            to_plot{k}=v_to_use(temp_idx);
            
        end
        
        if i==1
            %norm to WT
            v_temp=to_plot{1};
            for k=[1 3]
                to_plot{k}=to_plot{k}./v_temp;
            end

            v_temp=to_plot{2};
            for k=[2 4]
                to_plot{k}=to_plot{k}./v_temp;
            end
            
        else
            
            v_temp=to_plot{1};
            for k=1:4
                to_plot{k}=to_plot{k}./v_temp;
            end
            
        end
        
        for k=1:length(to_plot)
            
            v_mean(k)=mean(to_plot{k},'omitnan');
            v_sem(k)=std(to_plot{k},[],'omitnan')./sqrt(length(to_plot{k}));
            
        end
        
        if i==1
            plot_order=[1 3 2 4];
        else
            plot_order=1:4;
        end
        
        v_mean_all=[v_mean_all v_mean(plot_order)];
        v_sem_all=[v_sem_all v_sem(plot_order)];
        
        to_plot_all=[to_plot_all to_plot(plot_order)];
        
        
        if mod(j,2)==1
            v_merge=[];
            v_merge=[v_merge to_plot(plot_order)];
        elseif mod(j,2)==0
            v_merge=[v_merge to_plot(plot_order)];
        end
        
        m=m+1;
        
    end
    
end




v1=v_mean_all(17:32);
v2=v_sem_all(17:32);
v3=to_plot_all(17:32);

for i=condition_to_plot%1:2
    
    
    hold on
    
    to_plot1=v1(8*(i-1)+(1:4));
    to_plot2=v1(8*(i-1)+(5:8));
    to_plot3=v2(8*(i-1)+(1:4));
    to_plot4=v2(8*(i-1)+(5:8));
    
    to_plot5=v3(8*(i-1)+(1:4));
    to_plot6=v3(8*(i-1)+(5:8));
    
    bar([to_plot1; to_plot2]')
    errorbar((1:4)-0.15,to_plot1,to_plot3,'.k')
    errorbar((1:4)+0.15,to_plot2,to_plot4,'.k')
    
    for j=1:4
        [h p]=ttest2(to_plot5{j},to_plot6{j});
        text(j,1+0.05*j,num2str(p))
    end
    
    ylim([0.8 1.2])
    xlim([0.5 4.5])
    legend({'-rad','+rad'},'Location','Northwest')
    xticks(1:4)
    xtickangle(45)
    xticklabels(strain_names{2})
    ylabel('norm. to WT')
    
end










end



