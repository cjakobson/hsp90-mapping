function [] = plot_rictor_reconstruction(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



days={'20220623'};
times={'72h'};

input_data=readtable([dependency_directory days{1} 'radPhenotypingdata.csv']);
input_mat=table2array(input_data);
%image artifacts
input_mat(input_mat>3000)=nan;

n_plates=32;

conditions{3}={'glc-rad','glc+rad','rap-rad','rap+rad','glc+GdA','rap+GdA'};

strains{3}={'YJM975','RM11','YJM975 TSC11 1342Lys>Thr','RM11 TSC11 1342Thr>Lys'};


%rearrange to 96well

aIdxBase=1:2:24;
aIdx=aIdxBase;

bIdxBase=2:2:24;
bIdx=bIdxBase;

cIdxBase=25:2:48;
cIdx=cIdxBase;

dIdxBase=26:2:48;
dIdx=dIdxBase;


for j=1:7
    aIdx=[aIdx aIdxBase+48*j];
    bIdx=[bIdx bIdxBase+48*j];
    cIdx=[cIdx cIdxBase+48*j];
    dIdx=[dIdx dIdxBase+48*j];
end

for j=1:n_plates

    reorder_mat((1:96),j)=input_mat(aIdx,j);
    reorder_mat(((96+1):(2*96)),j)=input_mat(bIdx,j);
    reorder_mat(((2*96+1):(3*96)),j)=input_mat(cIdx,j);
    reorder_mat(((3*96+1):(4*96)),j)=input_mat(dIdx,j);

end


i=3;

plot_order=[1 3 2 4];

%plate offset for this strain/condition pair
r=1;
for j=1:length(conditions{i})
    
    plate_to_plot=6*(i-1)+j;

    temp_data=reorder_mat(:,plate_to_plot);

    clear v_to_plot
    for l=1:4

        v_to_plot{l}=temp_data(96*(l-1)+(1:96));
        
    end
    
    for l=1:4
        
        temp_data_mat(:,r)=v_to_plot{plot_order(l)};
        r=r+1;

    end
    
end


m=1;

%gather by genotype across conditions

%glc/YJM
idx_to_use{1}=[1 2 5 6 17 18];
%glc/RM
idx_to_use{2}=[3 4 7 8 19 20];
%rad/YJM
idx_to_use{3}=[9 10 13 14 21 22];
%rad/RM
idx_to_use{4}=[11 12 15 16 23 24];

for k=[3 4]%1:length(idx_to_use)
    clear v_to_plot v_mean v_sem
    for l=1:length(idx_to_use{k})

        v_to_plot{l}=temp_data_mat(:,idx_to_use{k}(l));
        v_mean(l)=mean(v_to_plot{l},'omitnan');
        v_sem(l)=std(v_to_plot{l},[],'omitnan')/sqrt(length(v_to_plot{l}));

    end

    for l=1:3
        v_sem(2*(l-1)+(1:2))=v_sem(2*(l-1)+(1:2))./v_mean(2*(l-1)+1);
        v_mean(2*(l-1)+(1:2))=v_mean(2*(l-1)+(1:2))./v_mean(2*(l-1)+1);
    end

    subplot(2,8,m)
    hold on
    bar(v_mean)
    errorbar(1:length(v_mean),v_mean,v_sem,'k.')
    %title([tempNames{k} times{n}])

    if k<=3
        ylim([0.7 1.1])
        
        for l=2:3

            [h,p]=ttest2(v_to_plot{2*(l-1)+2},v_to_plot{2});
            text(2*l,0.7+0.1*l,num2str(p))

        end
        
    else
        ylim([0.9 1.5])
        
        for l=2:3

            [h,p]=ttest2(v_to_plot{2*(l-1)+2},v_to_plot{2});
            text(2*l,1.1+0.1*l,num2str(p))

        end
        
    end



    

    m=m+1;



end
            




end





