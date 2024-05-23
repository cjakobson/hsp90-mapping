function [] = plot_vsrc(dependency_directory,output_directory)

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

days={'20220418'};
times={'48h'};

n=1;

inpute_data=readtable([dependency_directory days{n} 'vSrcPlatesdata.csv']);
input_mat=table2array(inpute_data);

n_plates=8;

conditions={'glc-rad','glc+rad','gal-rad','gal-rad'};

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

%just use first het instance
het_mat=reorder_mat((96+1):(2*96),:);

%just use first replicate
for i=1:length(conditions)
    to_plot{i}=het_mat(:,i);
end

hold on
easy_box(to_plot)
ylim([0 2500])
xticks(1:length(conditions))
xtickangle(45)
xticklabels(conditions)
ylabel('spot size')

[h p]=ttest2(to_plot{3},to_plot{4});
text(3.5,2000,num2str(p))
    


end

