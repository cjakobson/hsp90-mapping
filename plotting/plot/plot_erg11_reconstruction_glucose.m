function []= plot_erg11_reconstruction_glucose(condition_to_plot,plot_offset,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


days={'20220531'};
times={'72h'};


input_data=readtable([dependency_directory days{1} 'radPhenotypingdata.csv']);
input_mat=table2array(input_data);

n_plates=20;

conditions={'glc','mal','Ni','flc','teb'};
rad_conditions={'-rad','+rad'};

strains={'YJM975','RM11','RM11 ERG11 122014C>T','YJM975 ERG11 122014T>C'};


input_mat(input_mat<100)=nan;

%rearrange to 96well
%now to 96
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



%aggregate bars to reorder
vMean2=[];
vSem2=[];
vToPlot2=[];
for j=1:length(rad_conditions)

    for k=1%:length(strains)

        %put comparisons side by side
        plot_order=[1 4 2 3];

        plate_to_plot=4*(condition_to_plot-1)+2*(j-1)+k;

        temp_data=reorder_mat(:,plate_to_plot);

        for l=1:4

            v_to_plot{l}=temp_data(96*(l-1)+(1:96));
            v_mean(l)=mean(v_to_plot{l},'omitnan');
            v_sem(l)=std(v_to_plot{l},[],'omitnan')/sqrt(length(v_to_plot{l}));

        end

        %normalize to parental mean
        v_sem(plot_order(1:2))=v_sem(plot_order(1:2))./v_mean(plot_order(1));
        v_mean(plot_order(1:2))=v_mean(plot_order(1:2))./v_mean(plot_order(1));
        v_sem(plot_order(3:4))=v_sem(plot_order(3:4))./v_mean(plot_order(3));
        v_mean(plot_order(3:4))=v_mean(plot_order(3:4))./v_mean(plot_order(3));


        vMean2=[vMean2 v_mean];
        vSem2=[vSem2 v_sem];
        vToPlot2=[vToPlot2 v_to_plot];


    end

end


%separate by strain
idx_to_use=[1 4 5 8];
v_mean=vMean2(idx_to_use);
v_sem=vSem2(idx_to_use);
v_to_plot=vToPlot2(idx_to_use);


m=1;

subplot(2,8,plot_offset+m)
hold on
bar(v_mean)
errorbar(1:length(v_mean),v_mean,v_sem,'k.')
title([conditions{condition_to_plot} ' YJM975 ' times{1}])
xticks(1:length(strains{k}([1 4 1 4])))
xtickangle(45)
ylim([0.9 1.1])
plot(xlim,[1 1],':r')

m=m+1;





%RM11
idx_to_use=[2 3 6 7];
v_mean=vMean2(idx_to_use);
v_sem=vSem2(idx_to_use);
v_to_plot=vToPlot2(idx_to_use);

subplot(2,8,m+plot_offset)
hold on
bar(v_mean)
errorbar(1:length(v_mean),v_mean,v_sem,'k.')
title([conditions{condition_to_plot} ' RM11 ' times{1}])
xticks(1:length(strains{k}([1 4 1 4])))
xtickangle(45)
ylim([0.9 1.1])
plot(xlim,[1 1],':r')
ylabel('norm. to WT')


end



