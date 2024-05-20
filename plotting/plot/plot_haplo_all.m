function [] = plot_haplo_all(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

hap_table=readtable([dependency_directory 'pbio.2005130.s003.xlsx'],'Sheet','A');

essential_table=readtable([dependency_directory 'inviable_annotations_filtered_by_giaever.txt']);

mapping_input=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);
    
mapping_input=mapping_input(mapping_input.isQtn==1,:);



ess_idx=logical(ismember(mapping_input.gene1,essential_table.GeneSystematicName)+...
    ismember(mapping_input.gene2,essential_table.GeneSystematicName));

hsp_idx=logical(ismember(mapping_input.gene1,hap_table.ORF)+...
    ismember(mapping_input.gene2,hap_table.ORF));


v_delta_deltaZ=abs(mapping_input.deltaZbuffer-mapping_input.deltaZbaseline);

to_plot{1}=v_delta_deltaZ(logical(hsp_idx));
to_plot{2}=v_delta_deltaZ(logical(ess_idx.*~hsp_idx));
to_plot{3}=v_delta_deltaZ(logical(~ess_idx.*~hsp_idx));

hold on
easy_box(to_plot)
ylim([0 0.5])
for i=1:length(to_plot)
    text(i,0.45,num2str(length(to_plot{i})))
end
ylabel('\Delta\DeltaZ')
m=1;
for i=1:length(to_plot)
    for k=(i+1):length(to_plot)
        %[p h]=ranksum(toPlot{i},toPlot{k});
        %[p h]=ranksum(toPlot{i},toPlot{j});
        [h p]=ttest2(to_plot{i},to_plot{k});
        text((i+k)/2,0.4-0.025*m,num2str(p))
        m=m+1;
    end
end
tempLabels={'hap.','ess. not hap.','other'};
xticks(1:length(tempLabels))
xtickangle(45)
xticklabels(tempLabels)
title('all interactions')



end


