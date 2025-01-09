function [] = plot_hsp90_growth(plot_offset,dependency_directory,output_directory)

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
hsp82_idx=find(ismember(orf_names,'YPL240C'));

v_temp=input_mat(hsc82_idx,:);
v_temp=v_temp(strain_index>0);
v_hsc82(strain_index(strain_index>0))=v_temp;

v_temp=input_mat(hsp82_idx,:);
v_temp=v_temp(strain_index>0);
v_hsp82(strain_index(strain_index>0))=v_temp;




%load haploid growth
load([dependency_directory 'haploid_rad_trait.mat']);

v_no_rad=mean([trait{13} trait{14}],2);
v_rad=mean([trait{15} trait{16}],2);

subplot(2,4,plot_offset+1)
v1=v_hsc82';
v2=v_no_rad;
scatter(v1,v2,10,'k','filled')
xlim([1e5 2.5e5])
ylim([-3 3])
axis square
xlabel('Hsc82 abundance')
ylabel('growth in no stress')
[r p]=corr(v1,v2,'rows','complete');
text(1e5,-2,num2str(r))
text(1e5,-2.5,num2str(p))


subplot(2,4,plot_offset+2)
v1=v_hsc82';
v2=v_rad;
scatter(v1,v2,10,'k','filled')
xlim([1e5 2.5e5])
ylim([-3 3])
axis square
xlabel('Hsc82 abundance')
ylabel('growth in no stress+rad')
[r p]=corr(v1,v2,'rows','complete');
text(1e5,-2,num2str(r))
text(1e5,-2.5,num2str(p))


subplot(2,4,plot_offset+3)
v1=v_hsp82';
v2=v_no_rad;
scatter(v1,v2,10,'k','filled')
xlim([2e4 1e5])
ylim([-3 3])
axis square
xlabel('Hsp82 abundance')
ylabel('growth in no stress')
[r p]=corr(v1,v2,'rows','complete');
text(2e4,-2,num2str(r))
text(2e4,-2.5,num2str(p))


subplot(2,4,plot_offset+4)
v1=v_hsp82';
v2=v_rad;
scatter(v1,v2,10,'k','filled')
xlim([2e4 1e5])
ylim([-3 3])
axis square
xlabel('Hsp82 abundance')
ylabel('growth in no stress+rad')
[r p]=corr(v1,v2,'rows','complete');
text(2e4,-2,num2str(r))
text(2e4,-2.5,num2str(p))


%also against change in growth
subplot(2,4,plot_offset+5)
v1=v_hsc82';
v2=v_rad-v_no_rad;
scatter(v1,v2,10,'k','filled')
xlim([1e5 2.5e5])
ylim([-3 3])
axis square
xlabel('Hsc82 abundance')
ylabel('change in growth upon rad')
[r p]=corr(v1,v2,'rows','complete');
text(1e5,-2,num2str(r))
text(1e5,-2.5,num2str(p))


subplot(2,4,plot_offset+6)
v1=v_hsp82';
v2=v_rad-v_no_rad;
scatter(v1,v2,10,'k','filled')
xlim([2e4 1e5])
ylim([-3 3])
axis square
xlabel('Hsp82 abundance')
ylabel('change in growth upon rad')
[r p]=corr(v1,v2,'rows','complete');
text(2e4,-2,num2str(r))
text(2e4,-2.5,num2str(p))





end


