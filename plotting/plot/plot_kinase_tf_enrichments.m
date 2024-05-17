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

qThresh=0.05;


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




%TF targets from huge ChIPseq experiment
%https://www.nature.com/articles/s41586-021-03314-8#MOESM4

[num,txt]=xlsread([dependency_directory '41586_2021_3314_MOESM4_ESM.xlsx'],6);

vTF=txt(2,4:374);

vTFgenes=txt(9:5872,2);

vTFtype=txt(6,4:374);

clear txt

tfMat=num(9:5872,4:374);

clear num

toUse=ismember(vTFtype,'TF');

vTF=vTF(toUse);
tfMat=tfMat(:,toUse);


for i=1:length(vTF)
    
    tfTargets{i}=vTFgenes(logical(tfMat(:,i)));
    
end

for i=1:length(vTF)
    
    tfTargetOverlap{i,1}=intersect(tfTargets{i},qtn_genes);
    tfTargetOverlap{i,2}=intersect(tfTargets{i},qtl_genes);
    tfTargetOverlap{i,3}=intersect(tfTargets{i},all_segregating_genes);
    
    tempTable=table([length(tfTargetOverlap{i,1});length(qtn_genes)-length(tfTargetOverlap{i,1})],...
        [length(tfTargetOverlap{i,3});length(all_segregating_genes)-length(tfTargetOverlap{i,3})],...
        'VariableNames',{'interacting','all other'},'RowNames',{'client','not client'});
    [h,p,stats]=fishertest(tempTable);
    
    tfQtnPval(i)=p*length(vTF);
    
    if tfQtnPval(i)<0.05
        toOutput=table(tfTargetOverlap{i,1},'VariableNames',{'ORF'});
        writetable(toOutput,[vTF{i} '_qtnTargets.txt'])
    end
    
    
end


nMatTf=cellfun(@length,tfTargetOverlap);

for i=1:length(vTF)
        
    fMatTf(i,1)=nMatTf(i,1)/length(qtn_genes);
    fMatTf(i,2)=nMatTf(i,2)/length(qtl_genes);
    fMatTf(i,3)=nMatTf(i,3)/length(all_segregating_genes);  
     
end

tfQtnOddsRatio=fMatTf(:,1)./fMatTf(:,3);


qThresh=0.05;
toUse=tfQtnPval<qThresh;
qToUse=tfQtnPval(toUse);
tfToUse=vTF(toUse);



%sort by descending OR for merged fig
v1=tfQtnOddsRatio(toUse);
v2=qToUse;
v3=tfToUse';

[~,sortIdx]=sort(v1,'descend');

oddsForMerge=[oddsForMerge; v1(sortIdx)];
pForMerge=[pForMerge v2(sortIdx)];
labelsForMerge=[labelsForMerge; v3(sortIdx)];



hold on
bar(oddsForMerge)
ylim([0 3])
plot(xlim,[1 1],':r')
ylabel('QTN odds ratio')
xticks(1:length(oddsForMerge))
xtickangle(45)
xticklabels(labelsForMerge)
for i=1:length(pForMerge)
    if pForMerge(i)<qThresh
        h=text(i,2.2,num2str(pForMerge(i)));
        set(h,'Rotation',45);
    end
end




end

