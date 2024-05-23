function [] = plot_ase_all(dependency_directory,output_directory)

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

ase_input=readtable([dependency_directory 'radAseData.csv']);
mapping_input=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

%ase_idx=logical(ase_input.hasAseRad+ase_input.hasAseNoRad);
v1=ase_input.rmAfRad;%(ase_idx);
v2=ase_input.rmAfNoRad;%(ase_idx);

qtn_idx=mapping_input.isQtn==1;
v3=mapping_input.rmAfRad(qtn_idx);
v4=mapping_input.rmAfNoRad(qtn_idx);

hold on
scatter(v1./(1-v1),v2./(1-v2),10,'k','filled','MarkerFaceAlpha',0.5)
scatter(v3./(1-v3),v4./(1-v4),25,'r','filled','MarkerFaceAlpha',0.5)

axis square
set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([1e-2 1e2])
ylim(xlim)
plot(xlim,ylim,':r')
xlabel('RM allele ratio')
ylabel('RM allele ratio')



end