
clear

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',2)
figureCounter=1;


%filebase='/Users/cjakobson/';
filebase='/Users/christopherjakobson/';

code_directory=[filebase 'Documents/GitHub/hsp90-mapping/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/hsp90mapping/hsp90-mapping-dependencies/'];
%output_directory=[filebase 'Dropbox/JaroszLab/hsp90mapping/manuscript-plots/'];
output_directory=[filebase 'Dropbox/JaroszLab/hsp90mapping/revision-plots/'];

addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/parse'])
addpath([code_directory 'plotting/calculate'])
addpath([code_directory 'plotting/plot'])
addpath([code_directory 'data-prep'])



%load info on variants
variantInfo=readtable([dependency_directory 'variantInfoStructure.csv']);



%add whether Gong client, SGD genetic or physical interactor, etc
tempInput=readtable([dependency_directory 'GongTableS2.xlsx']);

chapLabels={'Hsp82','Hsc82','Hsp104','Hsp42','Ssa1','Ssa2','Ssa3','Ssa4',...
    'Sse1','Ydj1'};

clear upper
for i=1:length(chapLabels)
    
    tempData=table2array(tempInput(:,ismember(tempInput.Properties.VariableNames,chapLabels{i})));
    
    for j=1:length(tempData)
        
        if length(tempData{j})>0
            
            tempStr=strsplit(tempData{j},',');
            
            %make upper case to match bufferInput
            clientArray{i,j}=upper(tempStr{1});
            
        end
        
    end
    
end


for i=1:height(variantInfo)
    
    if ~isempty(variantInfo.common1{i})
        tempStr=strsplit(variantInfo.common1{i},';');
        tempCommon1{i}=tempStr{1};
    else
        tempCommon1{i}='NA';
    end
    if ~isempty(variantInfo.common2{i})
        tempStr=strsplit(variantInfo.common2{i},';');
        tempCommon2{i}=tempStr{1};
    else
        tempCommon2{i}='NA';
    end
    
end

tempCommon1(cellfun(@isempty,tempCommon1))={'NA'};
tempCommon2(cellfun(@isempty,tempCommon2))={'NA'};

for j=1:length(chapLabels)
        
    tempArray=clientArray(j,:);
    tempArray(cellfun(@isempty,tempArray))=[];

    v1=ismember(tempCommon1,tempArray);
    v2=ismember(tempCommon2,tempArray);

    vChap(:,j)=logical(v1+v2);
        
end


chapTable=array2table(vChap,'VariableNames',chapLabels);


toOutput=[variantInfo chapTable];
writetable(toOutput,[dependency_directory 'variantInfoChaperone.csv'])





