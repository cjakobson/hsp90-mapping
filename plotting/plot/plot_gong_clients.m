function [] = plot_gong_clients(dependency_directory,output_directory)


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

chap_sum_90_input=sum(logical(input_data.Hsp82+input_data.Hsc82));
chap_sum_70_input=sum(logical(input_data.Ssa1+input_data.Ssa2+...
    +input_data.Ssa3+input_data.Ssa4));
chap_sum_any_input=sum(logical(input_data.Hsp82+input_data.Hsc82+...
    input_data.Ssa1+input_data.Ssa2+...
    +input_data.Ssa3+input_data.Ssa4));

chap_frac_90_input=chap_sum_90_input./height(input_data)
chap_frac_70_input=chap_sum_70_input./height(input_data)
chap_frac_any_input=chap_sum_any_input./height(input_data)



background_data=readtable([dependency_directory 'variantInfoChaperone.csv']);

to_remove=unique(input_data.index);
background_data(to_remove,:)=[];

chap_sum_90_all=sum(logical(background_data.Hsp82+background_data.Hsc82));
chap_sum_70_all=sum(logical(background_data.Ssa1+background_data.Ssa2+...
    +background_data.Ssa3+background_data.Ssa4));
chap_sum_any_all=sum(logical(background_data.Hsp82+background_data.Hsc82+...
    background_data.Ssa1+background_data.Ssa2+...
    +background_data.Ssa3+background_data.Ssa4));

chap_frac_90_all=chap_sum_90_all./height(background_data)
chap_frac_70_all=chap_sum_70_all./height(background_data)
chap_frac_any_all=chap_sum_any_all./height(background_data)


%fishers
temp_table=table([chap_sum_any_input;height(input_data)-chap_sum_any_input],...
        [chap_sum_any_all;height(background_data)-chap_sum_any_all],...
        'VariableNames',{'interacting','all other'},'RowNames',{'client','not client'})
    [h,p,stats]=fishertest(temp_table)



end


