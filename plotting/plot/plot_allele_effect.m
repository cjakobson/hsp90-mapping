function [] = plot_allele_effect(locus_to_plot,condition_to_plot,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

%assume locus_to_plot doesn't have geometric offset

load([dependency_directory 'radFilename.mat'])
load([dependency_directory 'radTrait.mat'])

raw_size_data=readtable([dependency_directory 'radMeanStdData.csv']);


model_genotypes = parse_genotypes(dependency_directory,output_directory);
    
baseline_idx=ismember(filename,[condition_to_plot '-rad']);
buffer_idx=ismember(filename,[condition_to_plot '+rad']);

temp_baseline_trait=trait{baseline_idx};
temp_buffer_trait=trait{buffer_idx};

idx{1}=model_genotypes(:,locus_to_plot)==1;
idx{2}=model_genotypes(:,locus_to_plot)==-1;


v1=[mean(temp_baseline_trait(idx{1}),'omitnan')...
    mean(temp_baseline_trait(idx{2}),'omitnan')];
v2=[mean(temp_buffer_trait(idx{1}),'omitnan')...
    mean(temp_buffer_trait(idx{2}),'omitnan')];

v3=[std(temp_baseline_trait(idx{1}),[],'omitnan')/sqrt(sum(idx{1}))...
    std(temp_baseline_trait(idx{2}),[],'omitnan')/sqrt(sum(idx{2}))];
v4=[std(temp_buffer_trait(idx{1}),[],'omitnan')/sqrt(sum(idx{1}))...
    std(temp_buffer_trait(idx{2}),[],'omitnan')/sqrt(sum(idx{2}))];

%convert to raw size
temp_str1=[condition_to_plot '-rad'];
temp_str2=[condition_to_plot '+rad'];

mean1=raw_size_data.Var2(ismember(raw_size_data.conditions,temp_str1));
mean2=raw_size_data.Var2(ismember(raw_size_data.conditions,temp_str2));

std1=raw_size_data.Var3(ismember(raw_size_data.conditions,temp_str1));
std2=raw_size_data.Var3(ismember(raw_size_data.conditions,temp_str2));

v1=(v1*std1)+mean1;
v2=(v2*std2)+mean2;

v3=(v3*std1);
v4=(v4*std2);


hold on
bar([v1 v2])
errorbar([v1 v2],[v3 v4],'.k')
xlim([0.5 4.5])
ylim([0 Inf])
title([condition_to_plot ' ' num2str(locus_to_plot)])
%ylabel('norm. Z score')
ylabel('spot size')
n=1;
for k=1:2
    text(n,min(ylim),num2str(sum(idx{k})))
    n=n+1;
end
xticklabels({})




end