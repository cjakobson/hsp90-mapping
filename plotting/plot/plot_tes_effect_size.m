function [] =plot_tss_effect_size(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

input_table=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

input_table.deltaDeltaZ=abs(input_table.deltaZbuffer-input_table.deltaZbaseline);

table_to_use=parse_tss_tes(input_table,dependency_directory,output_directory);

qtn_idx=table_to_use.isQtn==1;


has_tag_idx=logical((table_to_use.hasAseRad==0).*(table_to_use.hasTagRad==1));

hold on
v1=table_to_use.tesDistGene1(logical(has_tag_idx.*qtn_idx.*(table_to_use.upstreamGene1==0)));
v2=abs(table_to_use.deltaDeltaZ(logical(has_tag_idx.*qtn_idx.*(table_to_use.upstreamGene1==0))));
scatter(v1,v2,50,'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor',grey)


ase_idx=table_to_use.hasAseRad==1;

v1=table_to_use.tesDistGene1(logical(ase_idx.*qtn_idx.*(table_to_use.upstreamGene1==0)));
v2=abs(table_to_use.deltaDeltaZ(logical(ase_idx.*qtn_idx.*(table_to_use.upstreamGene1==0))));
scatter(v1,v2,50,'filled','k')
axis square
xlim([-500 500])
ylim([0 0.75])
plot([0 0],ylim,':r')
ylabel('\Delta\DeltaZ')
xlabel('TES dist. (bp)')


end

