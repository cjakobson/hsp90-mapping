function [] = plot_mis_foldx(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

%just QTNs
input_data=input_data(input_data.isQtn==1,:);


[~,v_type]=variant_types(input_data.variantType);
v_type=v_type';


%v_regulatory=v_type==3;
v_mis=v_type==1;


to_plot{1}=abs(input_data.variantFoldX(v_mis));

%all other mis

background_data=readtable([dependency_directory 'variantInfoStructure.csv']);

[v_fraction_all,v_type_all]=variant_types(background_data.variantType);

to_plot{2}=abs(background_data.variantFoldX);

% (~ismember(1:height(background_data),...
%     input_data.index(v_mis))));



easy_box(to_plot)
%ylim([0 0.5])
% for i=1:length(to_plot)
%     text(i,0.4,num2str(length(to_plot{i})))
% end
% for i=2:length(to_plot)
%     [h p]=ttest2(to_plot{1},to_plot{i});
%     text(i-0.5,0.45,num2str(p))
% end
% ylabel('\DeltanTE')
% xticks(1:length(to_plot))
% xtickangle(45)
% xticklabels({'has ASE','no ASE'})%,'all segr. syn.'})
% title('synonymous variants')


end


