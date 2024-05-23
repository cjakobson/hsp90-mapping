function [] = plot_nfs1_reconstruction(conditions_to_plot,plot_offset,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



start_date=20230803;
start_time=13*60+48+24/60; 

plate_names={'a','b','c','d','e'};

n_scanners=2;

condition_names={'glc-rad','dox-rad','flc-rad','rap-rad','teb-rad',...
    'glc+rad','dox+rad','flc+rad','rap+rad','teb+rad'};

strain_names={'6635_YJM975','6649_RM11','8525_RM11_NCP1-1','9255_YJM975_NFS1-1'};




%get filenames
to_read=dir([dependency_directory 'nfs1rad/gitter']);

%get sga data
m=1;
for i=1:length(to_read)
    
    temp_name=to_read(i).name;

    if length(temp_name)>3

        if strcmp(temp_name((end-2):end),'dat')

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
            sga_mat{m}=readtable([dependency_directory 'nfs1rad/gitter/' temp_name]);

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



v_merge=[];

m=1;
for j=conditions_to_plot
        
    temp_mat=reorder_mat{j};

    fwhm_idx=72;

    v_to_use=temp_mat(:,fwhm_idx);

    clear to_plot
    for k=1:length(strain_names)

        temp_idx=384*(k-1)+(1:384);
        to_plot{k}=v_to_use(temp_idx);

    end

    v_merge=[v_merge to_plot];

    %norm to WT
    v_temp=to_plot{1};
    for k=[1 4]
        to_plot{k}=to_plot{k}./v_temp;
    end

    v_temp=to_plot{2};
    for k=[2 3]
        to_plot{k}=to_plot{k}./v_temp;
    end


    for k=1:length(to_plot)

        v_mean(k)=mean(to_plot{k},'omitnan');
        v_sem(k)=std(to_plot{k},[],'omitnan')./sqrt(length(to_plot{k}));

    end

    strain_to_use=[1 4];

    subplot(2,16,plot_offset+m)
    hold on
    bar(v_mean(strain_to_use))
    errorbar(v_mean(strain_to_use),v_sem(strain_to_use),'.k')
    title(condition_names{j})
    xticks(1:length(strain_names(strain_to_use)))
    xtickangle(45)
    xticklabels(strain_names(strain_to_use))
    ylim([0.6 1.4])
    plot(xlim,[1 1],':k')
    ylabel('norm. to WT')

    [h p]=ttest2(to_plot{strain_to_use(1)},to_plot{strain_to_use(2)});
    text(2,1.3,num2str(p))


    m=m+1;

end



end



