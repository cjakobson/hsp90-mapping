function [] = plot_erg11_crispri_norm_glc(strain_to_plot,condition_to_plot,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


exp_name='20230929sigmoid';



condition_names={'glc+aTc0','glc+aTc62','glc+aTc125','glc+aTc250','glc+aTc500','glc+aTc1000'...
    'glc+aTc0+rad','glc+aTc62+rad','glc+aTc125+rad','glc+aTc250+rad','glc+aTc500+rad','glc+aTc1000+rad',...
    'flc25+aTc0','flc25+aTc62','flc25+aTc125','flc25+aTc250','flc25+aTc500','flc25+aTc1000'...
    'flc25+aTc0+rad','flc25+aTc62+rad','flc25+aTc125+rad','flc25+aTc250+rad','flc25+aTc500+rad','flc25+aTc1000+rad',...
    'flc100+aTc0','flc100+aTc62','flc100+aTc125','flc100+aTc250','flc100+aTc500','flc100+aTc1000'...
    'flc100+aTc0+rad','flc100+aTc62+rad','flc100+aTc125+rad','flc100+aTc250+rad','flc100+aTc500+rad','flc100+aTc1000+rad'};
    

strain_names={'YJM975 no guide','RM11 no guide','YJM975 antiERG11','RM11 antiERG11'};

v_time={'48h'};


n_images=12;
n_plates=36;

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

        if exist([dependency_directory 'sigmoid-glc/gitter/' temp_name])

            sga_mat{m}=readtable([dependency_directory 'sigmoid-glc/gitter/' temp_name]);

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
        
        column_offset=(i-1)*length(condition_names);
        idx_to_use=column_offset+j;
        
        clear to_plot
        
        for k=1:length(strain_names)
            
            to_plot{k}=reorder_mat(384*(k-1)+(1:384),idx_to_use);
            
        end
        
        
        %normalize to no edit
        v_temp=to_plot{1};
        for ii=[1 3]%1:length(toPlot)
            to_plot{ii}=to_plot{ii}./v_temp;
        end
        v_temp=to_plot{2};
        for ii=[2 4]%1:length(toPlot)
            to_plot{ii}=to_plot{ii}./v_temp;
        end

        
        for ii=1:length(to_plot)
            
            to_plot{ii}(isinf(to_plot{ii}))=nan;
            
            v_mean(ii)=mean(to_plot{ii},'omitnan');
            v_sem(ii)=std(to_plot{ii},[],'omitnan')./sqrt(length(to_plot{ii}));
        
        end
        
        
        v_mean_all=[v_mean_all v_mean];
        v_sem_all=[v_sem_all v_sem];
        to_plot_all=[to_plot_all to_plot];
        

        
    end
    
end



%separate out RM and YJM
temp_labels1={'YJM975 antiERG11','RM11 antiERG11'};
temp_labels2={'glucose','fluconazole25','fluconazole100'};
temp_labels3={'0ng/uL','62ng/uL','125ng/uL','250ng/uL','500ng/uL','1000ng/uL'};

m=1;

offset=6;

for i=strain_to_plot%1:length(temp_labels1)
    
    v_mean_to_plot=v_mean_all((2+i):4:length(v_mean_all));
    v_sem_to_plot=v_sem_all((2+i):4:length(v_sem_all));
    v_to_plot_all=to_plot_all((2+i):4:length(v_sem_all));
    
    for j=condition_to_plot%1:length(temp_labels2)
        
        temp_idx=2*offset*(j-1)+(1:(2*offset));
        
        v1=v_mean_to_plot(temp_idx);
        v2=v_sem_to_plot(temp_idx);
        v3=v_to_plot_all(temp_idx);
        
        hold on
        %-rad
        bar([v1(1:offset); v1((offset+1):(2*offset))]')
        errorbar((1:offset)-0.15,v1(1:offset),v2(1:offset),'.k')
        %+rad
        errorbar((1:offset)+0.15,v1((1+offset):(2*offset)),v2((1+offset):(2*offset)),'.k')
        
        %p vals for rad vs not
        for ii=1:offset
            [h p]=ttest2(v3{ii},v3{ii+offset});
            text(ii,1+0.01*ii,num2str(p))
        end
        
        title([temp_labels1{i} ' ' temp_labels2{j}])
        legend({'-rad','+rad'})
        xticks(1:length(temp_labels3))
        xtickangle(45)
        xticklabels(temp_labels3)
        xlim([1.5 6.5])
        %ylim([0.4 1.2])
        ylabel('growth relative to non-targeting guide')
        xlabel('aTc conc. (CRISPRi expression)')
        axis square
        
        m=m+1;
            
    end
    
end



end


