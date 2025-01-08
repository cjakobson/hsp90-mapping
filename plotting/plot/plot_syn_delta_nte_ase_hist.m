function [] = plot_syn_delta_nte_ase_hist(dependency_directory,output_directory)


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

v_buffer=abs(input_data.deltaZbuffer)>=abs(input_data.deltaZbaseline);

%v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);
%v_delta_cai=abs(input_data.variantDeltaCai);
v_delta_nte=abs(input_data.refnTE-input_data.altnTE);

v_has_ase=logical(input_data.hasAseRad);

v_has_tag=logical(input_data.hasTagNoRad+input_data.hasTagRad);

[~,v_type]=variant_types(input_data.variantType);
v_type=v_type';


%v_regulatory=v_type==3;
v_regulatory=v_type==2;


to_plot{1}=v_delta_nte(logical(v_has_ase.*v_regulatory));
to_plot{2}=v_delta_nte(logical(v_has_tag.*~v_has_ase.*v_regulatory));

%all other syn

background_data=readtable([dependency_directory 'variantInfoStructure.csv']);

[v_fraction_all,v_type_all]=variant_types(background_data.variantType);

v_nte_all=abs(background_data.refnTE-background_data.altnTE);
%to_plot{3}=v_nte_all(v_type_all==2);


hold on
histogram(to_plot{1},0:0.1:1,'Normalization','probability')
histogram(to_plot{2},0:0.1:1,'Normalization','probability')
axis square
legend({'has ASE','no ASE'})%,'all segr. syn.'})
title('synonymous variants')
xlabel('\DeltanTE')
ylim([0 1])
ylabel('relative frequency')



end


