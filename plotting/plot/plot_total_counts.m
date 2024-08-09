function []= plot_total_counts(locus_to_plot,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


ase_data=readtable([dependency_directory 'radAseData.csv']);

%ref no rad
to_plot{1}=ase_data{locus_to_plot,16:18}+ase_data{locus_to_plot,22:24};
%alt no rad
%to_plot{2}=ase_data{locus_to_plot,22:24};


%ref rad
to_plot{2}=ase_data{locus_to_plot,13:15}+ase_data{locus_to_plot,19:21};
%alt rad
%to_plot{4}=ase_data{locus_to_plot,19:21};


for i=1:length(to_plot)
    v_mean(i)=mean(to_plot{i},'omitnan');
end

hold on
bar(v_mean)
for i=1:length(to_plot)
    scatter(i*ones(1,length(to_plot{i})),to_plot{i},10,'k','filled')
end
xlim([0.5 2.5])
ylabel('est. counts')


end

