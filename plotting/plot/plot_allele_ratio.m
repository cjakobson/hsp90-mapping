function []= plot_allele_ratio(locus_to_plot,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


ase_data=readtable([dependency_directory 'radAseData.csv']);

%ref no rad
to_plot{1}=ase_data{locus_to_plot,16:18};
%alt no rad
to_plot{2}=ase_data{locus_to_plot,22:24};


%ref rad
to_plot{3}=ase_data{locus_to_plot,13:15};
%alt rad
to_plot{4}=ase_data{locus_to_plot,19:21};

%rearrange if RM is alt
if ase_data.altAllele(locus_to_plot)==ase_data.rmAllele(locus_to_plot)
    v_temp=to_plot;
    to_plot{1}=v_temp{2};
    to_plot{2}=v_temp{1};
    
    to_plot{3}=v_temp{4};
    to_plot{4}=v_temp{3};
end

for i=1:2
    ratio_to_plot{i}=to_plot{2*(i-1)+1}./to_plot{2*(i-1)+2};
    v_ratio(i)=mean(ratio_to_plot{i},'omitnan');
end

hold on
bar(v_ratio)
for i=1:length(ratio_to_plot)
    scatter(i*ones(1,length(ratio_to_plot{i})),ratio_to_plot{i},10,'k','filled')
end
xlim([0 3])
set(gca,'YScale','log')
ylim([1/3 3])
plot(xlim,[1 1],':r')



end