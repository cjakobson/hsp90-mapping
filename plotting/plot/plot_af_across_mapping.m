function [] = plot_af_across_mapping(plot_offset,dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

variant_info=readtable([dependency_directory 'variantInfoStructure.csv']);

hsp90_input=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);
no_rad_input=readtable([dependency_directory 'linear_no_rad_fdr_0.05.csv']);
rad_input=readtable([dependency_directory 'linear_rad_fdr_0.05.csv']);

hsp90_qtn_idx=hsp90_input.isQtn==1;
v_buffer=abs(hsp90_input.deltaZbuffer)>=abs(hsp90_input.deltaZbaseline);
v_potentiate=abs(hsp90_input.deltaZbuffer)<abs(hsp90_input.deltaZbaseline);

clear to_plot
%to_plot{1}=abs(hsp90_input.deltaZbaseline(logical(hsp90_qtn_idx.*v_buffer)));
to_plot{1}=hsp90_input.nMinor(logical(hsp90_qtn_idx.*v_buffer));
%to_plot{2}=abs(hsp90_input.deltaZbaseline(logical(hsp90_qtn_idx.*v_potentiate)));
to_plot{2}=hsp90_input.nMinor(logical(hsp90_qtn_idx.*v_potentiate));

%to_plot{3}=abs(no_rad_input.deltaZbaseline(no_rad_input.isQtn==1));
to_plot{3}=no_rad_input.nMinor(no_rad_input.isQtn==1);

%to_plot{4}=abs(rad_input.deltaZbaseline(rad_input.isQtn==1));
%to_plot{4}=rad_input.nMinor(rad_input.isQtn==1);

to_plot{4}=variant_info.nMinor;


subplot(2,8,plot_offset+1)
hold on
easy_box(to_plot)
%ylim([-0.1 1])
temp_labels={'buffered','potentiated','found in no rad','all segregating'};
xticks(1:length(temp_labels))
xticklabels(temp_labels)
xtickangle(45)
title('AF in 1K')
ylabel('strains with minor allele')
for i=1:length(to_plot)
    text(i,0,num2str(length(to_plot{i})))
end
for i=2:length(to_plot)
    [h p]=ttest2(to_plot{1},to_plot{i});
    text((i+1)/2,350+25*i,num2str(p))
    plot([1 i],[350+25*i 350+25*i],'-k')
end
for i=3:length(to_plot)
    [h p]=ttest2(to_plot{2},to_plot{i});
    text((i+1)/2,250+25*i,num2str(p))
    plot([2 i],[250+25*i 250+25*i],'-k')
end
xlim([0.5 length(to_plot)+0.5])


end


