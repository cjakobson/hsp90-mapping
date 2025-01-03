function [] = plot_hub_correlation(common_input,systematic_input,dependency_directory,output_directory)

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


[input_mat,strain_names,ydj_names,strain_merge_idx,rm_idx,yjm_idx,f6_idx,sgrp_idx,...
    orf_names,strain_index]=parse_raw_abundance(dependency_directory,output_directory);

hsc82_idx=find(ismember(orf_names,'YMR186W'));
query_idx=find(ismember(orf_names,systematic_input));

v1=input_mat(hsc82_idx,f6_idx);
v2=input_mat(query_idx,f6_idx);

scatter(v1,v2,10,'k','filled')
axis square
ylabel([common_input ' abundance'])
xlabel('Hsc82 abundance')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% ylabel('# of interactors')
% xlabel('est. protein abundance')
[r p]=corr(v1',v2','rows','complete');
temp_ylim=ylim;
text(2e5,temp_ylim(2),['r = ' num2str(r)])
text(2e5,0.95*temp_ylim(2),['p = ' num2str(p)])




end


