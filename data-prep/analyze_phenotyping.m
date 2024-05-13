function []=analyze_phenotyping(dependency_directory,output_directory)


figureCounter=1; %figure counter

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)


expName={'20211214phenotyping','20211215phenotyping','20211216phenotyping','20211217phenotyping'};
expTime={'24h','48h','72h','96h'};

nImages=80;
nPlates=240;

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


platesPerCondition=10;
m=1;
for i=1:length(conditions)
    
    conditionMat{i}=[];
    
    for j=1:platesPerCondition
        conditionMat{i}=[conditionMat{i}; fullMat(:,m)];
        m=m+1;
    end 
    
    %remove outlying plates (gitter errors)
    tempMean=mean(conditionMat{i},'omitnan');
    tempStd=std(conditionMat{i},[],'omitnan');
    for j=1:platesPerCondition
        tempMat=conditionMat{i}((1536*(j-1)+1):(1536*j));
        tempZ=mean((tempMat-tempMean)./tempStd,'omitnan');
        if tempZ>1
            conditionMat{i}((1536*(j-1)+1):(1536*j))=nan;
        end
    end
    
end




%z score by conditions
for i=1:length(conditions)
    
    vTemp=conditionMat{i};
    zMat(:,i)=(vTemp-mean(vTemp,'omitnan'))./std(vTemp,[],'omitnan');
    
end

zMat(abs(zMat)>5)=nan;


%output trait and filename for mapping
filename=conditions;
for i=1:length(filename)
    trait{i}=zMat(:,i);
end

%plate B+1 on 2021/11/15 had gitter issue--remove
trait{28}(1:nSpots)=nan;

%also deltaZ for buffer mapping
m=length(trait)+1;
for i=1:2:length(filename)
    v1=trait{i};
    v2=trait{i+1};
    trait{m}=v1-v2;
    
    tempStr=strsplit(filename{i},'-');
    filename{m}=[tempStr{1} '_delta'];
    
    m=m+1;
end

save([output_directory 'radFilename.mat'],'filename')
save([output_directory 'radTrait.mat'],'trait')

to_output=table(filename);
writetable(to_output,[output_directory 'sample_names.txt'])

%make some plots with raw colony size data
m=1;
for i=1:2:length(conditions)
    
    meanNoRad(i)=mean(conditionMat{i},'omitnan');
    stdNoRad(i)=std(conditionMat{i},[],'omitnan');
    meanToOutput(m)=meanNoRad(i);
    stdToOutput(m)=stdNoRad(i);
    m=m+1;
    
    meanRad(i)=mean(conditionMat{i+1},'omitnan');
    stdRad(i)=std(conditionMat{i+1},[],'omitnan');
    meanToOutput(m)=meanRad(i);
    stdToOutput(m)=stdRad(i);
    m=m+1;
    
end

%output these data to add %growth changes to mapping output
toOutput=table(conditions,meanToOutput',stdToOutput');
writetable(toOutput,[output_directory 'radMeanStdData.csv'])



end





