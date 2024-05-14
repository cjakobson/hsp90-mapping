function [] = plot_1K_effect(locus_to_plot,condition_to_plot,time_to_plot,dependency_directory,output_directory)



set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;





expName={'20221030phenotyping','20221031phenotyping','20221101phenotyping','20221102phenotyping'};
expTime={'24h','48h','72h','96h'};

nImages=34;
nPlates=102;

plateNames={'a','b','c'};

inputData=[];
for i=1:length(expName)
    tempTable=readtable([dependency_directory expName{i} 'data.csv']);
    inputData=[inputData table2array(tempTable)];
end

[num txt]=xlsread([dependency_directory '20211209 radicicol rescreen plate key_FINAL.xlsx']);



tempConditions1=txt(2:13,2);
tempConditions2=num(1:12,1);
tempConditions3=txt(2:13,4);

for i=1:length(tempConditions1)
    tempConditions{i}=[tempConditions1{i} '_' num2str(tempConditions2(i))...
        tempConditions3{i}];
end

m=1;
for i=1:length(expTime)
    for j=1:length(tempConditions)
        radString={'-rad','+rad'};
        for k=1:2
            conditions{m}=[expTime{i} ' ' tempConditions{j} radString{k}];
            m=m+1;
        end
    end
end
conditions=conditions';


nSpots=1536;

fullMat=[];
for i=1:(nPlates*length(expTime))
        
    
    tempMat{i}=inputData(:,i);
    
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

    reorderMat{i}=nan(size(tempMat{i}));

    reorderMat{i}((1:384))=tempMat{i}(a1idx);
    reorderMat{i}(((384+1):(2*384)))=tempMat{i}(a2idx);
    reorderMat{i}(((2*384+1):(3*384)))=tempMat{i}(b1idx);
    reorderMat{i}(((3*384+1):(4*384)))=tempMat{i}(b2idx);

    fullMat=[fullMat reorderMat{i}];
    
end

fullMat(fullMat<25)=nan;
fullMat(fullMat>2000)=nan;


%trim to just 1K
plateOffset=102;
plateRange=1:48;

toKeep=[];
for i=1:length(expTime)
    toKeep=[toKeep (i-1)*plateOffset+plateRange];
end

fullMat=fullMat(:,toKeep);

nPlates=48;

%calculate change in growth for each replicate
for i=2:2:(nPlates*length(expTime))
    
    deltaMat(:,i/2)=(fullMat(:,i)-fullMat(:,i-1))./fullMat(:,i-1);
    
end

for i=2:2:nPlates*length(expTime)/2
    
    meanMat(:,i/2)=mean(deltaMat(:,[i-1 i]),2,'omitnan');
    
end

%get strain info
strainInfo=readtable([dependency_directory '1011 Genomes_Sace_strains_matrix_positions_384.xlsx']);
m=1;
for i=1:3   %plates
    
    tempIdx1=ismember(strainInfo.Matrix_384,['M' num2str(i)]);
    
    for j=1:16  %rows
        
        tempIdx2=ismember(strainInfo.row_384,j);
        
        for k=1:24  %columns
            
            tempIdx3=ismember(strainInfo.col_384,k);
            
            tempIdx=find(tempIdx1.*tempIdx2.*tempIdx3);
            
            if length(tempIdx)>0
                
                name1{m}=strainInfo.Strain_Name_4SRA{tempIdx};
                name2{m}=strainInfo.Standardized_name{tempIdx};
                
            else
                
                name1{m}='NA';
                name2{m}='NA';
                
            end
            
            m=m+1;
            
        end
        
    end
        
end


%should recast this to be deltaGrowth -/+ rad
%calculate change in growth for each replicate
for i=1:(nPlates*length(expTime)/4)
    
    idx1=(i-1)*4+1;
    idx2=(i-1)*4+3;
    
    noRadMeanMat(:,i)=mean([fullMat(:,idx1) fullMat(:,idx2)],2);
    
    idx1=(i-1)*4+2;
    idx2=(i-1)*4+4;
    
    radMeanMat(:,i)=mean([fullMat(:,idx1) fullMat(:,idx2)],2);
    
end


   

condition_to_use=find(ismember(tempConditions1,condition_to_plot));
time_to_use=find(ismember(expTime,time_to_plot));

%minor allele frequency from 1k genomes
load([dependency_directory '1002data.mat'])

tempHetMat=minGenotype~=minGenotype2;

minGenotype(tempHetMat)=-1;

for i=1:length(locus_to_plot)

    altIdx=logical((minGenotype(locus_to_plot(i),:)==1)+(minGenotype(locus_to_plot(i),:)==-1));
    refIdx=logical((minGenotype(locus_to_plot(i),:)==0));
    
    altStrains{i}=strainString(altIdx);
    refStrains{i}=strainString(refIdx);

end


%normalize to plot for manuscript
for j=time_to_use%1:length(expTime)
    
    columnOffset=nPlates/4*(j-1);
    
    for i=condition_to_use%1:length(tempConditions1)

        conditionIdx=i;%conditionToUse(i);
        hold on
        toPlot{2}=noRadMeanMat(1:1152,conditionIdx+columnOffset);

        altIdx=find(ismember(name2,altStrains{1}));

        toPlot{1}=toPlot{2}(altIdx);
        toPlot{2}(altIdx)=[];
        
        toPlot{4}=radMeanMat(1:1152,conditionIdx+columnOffset);
        
        toPlot{3}=toPlot{4}(altIdx);
        toPlot{4}(altIdx)=[];
        
        clear vMean vSem
        for k=1:length(toPlot)
            vMean(k)=median(toPlot{k},'omitnan');
            vSem(k)=std(toPlot{k},[],'omitnan')./sqrt(length(toPlot{k}));
        end
        
        vSem(1:2)=vSem(1:2)./vMean(1);
        vMean(1:2)=vMean(1:2)./vMean(1);
        
        vSem(3:4)=vSem(3:4)./vMean(3);
        vMean(3:4)=vMean(3:4)./vMean(3);
            
        
        bar(vMean)
        errorbar(1:length(vMean),vMean,vSem,'k.')

        ylabel('relative growth')
        title([condition_to_plot ' ' time_to_plot])
        tempLabels={'ref-rad','alt-rad','ref+rad','alt+rad'};
        xticks(1:length(tempLabels))
        xticklabels(tempLabels)
        xtickangle(45)
        [h p]=ttest2(toPlot{2},toPlot{4});
        text(3,1.7,num2str(p))
        xlim([0 5])

        

    end
    

end




end



