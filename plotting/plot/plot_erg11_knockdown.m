function [] = plot_erg11_knockdown(dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;


fl_input=readtable([dependency_directory '20250225erg11_fusion_knockdown.txt']);

fl_data=fl_input([10:2:24 56:2:70],1);

for i=1:16

    temp_str=strsplit(fl_data.ReadHeight{i},'\t');

    for j=1:12

        mng_mat(i,j)=str2num(temp_str{j+1});

    end

end

%rearrange plates
reorder_mat=[mng_mat(1:8,1:12) mng_mat(9:16,1:12)];

strain_names={'wt scramble','wt erg11',...
    'erg11::mNG scramble','erg11::mNG erg11'};

concs={'1000','500','250','125','62'};%,'31','16','8'};


for i=1:length(concs)

    v_temp=reorder_mat(i,:);

    for j=1:length(strain_names)

        temp_idx=j:4:length(v_temp);

        data_array{i,j}=v_temp(temp_idx);

    end

    %norm KD to scramble for fusion
    norm_array{i}=data_array{i,4}./data_array{i,3};
    v_mean(i)=mean(norm_array{i});

end


hold on
plot(1:length(v_mean),v_mean)
for i=1:length(norm_array)
    scatter(i,norm_array{i},25,'k','filled')
end
axis square
xlabel('[aTc]')
ylabel('rel. fluorescence')
xticks(1:length(concs))
xticklabels(concs)
ylim([0.5 1])
title('ERG11::mNeonGreen level')
xlim([0.5 length(concs)+0.5])

end


