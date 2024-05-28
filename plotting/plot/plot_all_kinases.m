function [] = plot_kinase_tf_enrichments(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

qtn_idx=input_data.isQtn==1;
temp_genes=unique([input_data.gene1(qtn_idx);input_data.gene2(qtn_idx)]);
temp_genes(cellfun(@isempty,temp_genes))=[];

qtn_genes=temp_genes;


temp_genes=unique([input_data.gene1;input_data.gene2]);
temp_genes(cellfun(@isempty,temp_genes))=[];

qtl_genes=temp_genes;


variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

temp_genes=unique([variant_info.gene1;variant_info.gene2]);
temp_genes(cellfun(@isempty,temp_genes))=[];

all_segregating_genes=temp_genes;
all_segregating_genes(ismember(all_segregating_genes,qtn_genes))=[];




v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);


%kinase data from MOTIPS/https://www.science.org/doi/10.1126/scisignal.2000482
%(dataset 2)

tempNames=dir([dependency_directory 'Dataset S2']);

m=1;
for i=4:length(tempNames)
    
    kinases{m}=tempNames(i).name;
    m=m+1;
    
end

kinases=kinases';

%threshold on likelihood to include
thresh=1;
for i=1:length(kinases)
    
    tempNames=dir([dependency_directory 'Dataset S2/' kinases{i}]);
    toGet=tempNames(3).name;
    
    tempData=readtable([dependency_directory 'Dataset S2/' kinases{i} '/' toGet],'FileType','text');
    
    targets{i}=tempData.ORF(tempData.Likelihood>=thresh);
    
    targetOverlap{i,1}=intersect(targets{i},qtn_genes);
    targetOverlap{i,2}=intersect(targets{i},qtl_genes);
    targetOverlap{i,3}=intersect(targets{i},all_segregating_genes);
    
    tempTable=table([length(targetOverlap{i,1});length(qtn_genes)-length(targetOverlap{i,1})],...
        [length(targetOverlap{i,3});length(all_segregating_genes)-length(targetOverlap{i,3})],...
        'VariableNames',{'interacting','all other'},'RowNames',{'client','not client'});
    [h,p,stats]=fishertest(tempTable);
    
    qtnPval(i)=p*length(kinases);
    
    
end

nMat=cellfun(@length,targetOverlap);

for i=1:length(kinases)
        
    fMat(i,1)=nMat(i,1)/length(qtn_genes);
    fMat(i,2)=nMat(i,2)/length(qtl_genes);
    fMat(i,3)=nMat(i,3)/length(all_segregating_genes);  
     
end

qtnOddsRatio=fMat(:,1)./fMat(:,3);

%qThresh=0.05;
qThresh=Inf;


toUse=qtnPval<qThresh;
qToUse=qtnPval(toUse);
kinaseToUse=kinases(toUse);


v1=qtnOddsRatio(toUse);
v2=qToUse;
v3=kinaseToUse;

[~,sortIdx]=sort(v1,'descend');

oddsForMerge=[];
pForMerge=[];
labelsForMerge=[];


oddsForMerge=[oddsForMerge; v1(sortIdx)];
pForMerge=[pForMerge v2(sortIdx)];
labelsForMerge=[labelsForMerge; v3(sortIdx)];



hold on
bar(oddsForMerge)
%ylim([0 3])
set(gca,'YScale','log')
ylim([0.5 3])
plot(xlim,[1 1],':r')
ylabel('QTN odds ratio')
xticks(1:length(oddsForMerge))
xtickangle(45)
xticklabels(labelsForMerge)
for i=1:length(pForMerge)
    if pForMerge(i)<0.05
        h=text(i,2.2,num2str(pForMerge(i)));
        set(h,'Rotation',45);
    end
end




end

