function [] = plot_effect_size_line_crossing(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)



blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



input_data=readtable([dependency_directory 'linear_hsp90_fdr_0.05.csv']);

v_buffer=abs(input_data.deltaZbuffer)>=abs(input_data.deltaZbaseline);

v_line_crossing=(input_data.deltaZbuffer.*input_data.deltaZbaseline)<0;

v_delta_delta_z=abs(input_data.deltaZbuffer-input_data.deltaZbaseline);

to_plot{1}=v_delta_delta_z(logical(v_buffer.*~v_line_crossing));
to_plot{2}=v_delta_delta_z(logical(v_buffer.*v_line_crossing));

to_plot{3}=v_delta_delta_z(logical(~v_buffer.*~v_line_crossing));
to_plot{4}=v_delta_delta_z(logical(~v_buffer.*v_line_crossing));


hold on
easy_box(to_plot)
xlim([0.5 4.5])
[h p]=ttest2(to_plot{1},to_plot{2});
text(1.5,0.4,num2str(p))
[h p]=ttest2(to_plot{3},to_plot{4});
text(3.5,0.4,num2str(p))
set(gca,'YScale','log')
ylim([0.01 1])
ylabel('\Delta\DeltaZ')
for i=1:length(to_plot)
    text(i,0.9,num2str(length(to_plot{i})))
end


end


