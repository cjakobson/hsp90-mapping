function [] = plot_variant_types(dependency_directory,output_directory)



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

[v_fraction_hsp90,v_type_hsp90]=variant_types(input_data.variantType);


%contributions by type
for i=1:3

    temp_idx=v_type_hsp90==i;
    v_contribution(i)=sum(input_data.varExp(temp_idx));

end

v_contribution=v_contribution./sum(v_contribution)

background_data=readtable([dependency_directory 'variantInfoStructure.csv']);

[v_fraction_all,v_type_all]=variant_types(background_data.variantType);



variant_labels={'protein-altering','synonymous','outside ORFs'};

for i=1:length(variant_labels)

    to_plot_hsp90(i)=sum(v_type_hsp90==i);
    to_plot_all(i)=sum(v_type_all==i);

end

to_plot_hsp90
to_plot_all
to_plot_hsp90./to_plot_all

bar([v_fraction_hsp90;v_fraction_all]')
axis square
xtickangle(45)
xticklabels(variant_labels)
ylabel('relative frequency')
legend({'Hsp90-dep.','all segr.'})
xlim([0.5 length(variant_labels)+0.5])



end

