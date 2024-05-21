function [] = plot_ase_reproducibility(dependency_directory,output_directory)

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;

ase_data=readtable([dependency_directory 'radAseData.csv']);

%ref no rad
to_plot{1}=ase_data{:,16:18};
%alt no rad
to_plot{2}=ase_data{:,22:24};


%ref rad
to_plot{3}=ase_data{:,13:15};
%alt rad
to_plot{4}=ase_data{:,19:21};

%rearrange if RM is alt
v_temp=to_plot;

to_plot{1}(ase_data.altAllele==ase_data.rmAllele,:)=...
    v_temp{2}(ase_data.altAllele==ase_data.rmAllele,:);
to_plot{2}(ase_data.altAllele==ase_data.rmAllele,:)=...
    v_temp{1}(ase_data.altAllele==ase_data.rmAllele,:);

to_plot{3}(ase_data.altAllele==ase_data.rmAllele,:)=...
    v_temp{4}(ase_data.altAllele==ase_data.rmAllele,:);
to_plot{4}(ase_data.altAllele==ase_data.rmAllele,:)=...
    v_temp{3}(ase_data.altAllele==ase_data.rmAllele,:);



m=1;
for i=1:2
    
    fraction_to_plot{i}=to_plot{2*(i-1)+1}./(to_plot{2*(i-1)+1}+to_plot{2*(i-1)+2});
    
    for j=1:3
        
        for k=(j+1):3
            
            v1=fraction_to_plot{i}(:,j);
            v2=fraction_to_plot{i}(:,k);
            
            subplot(2,3,m)
            scatter(v1,v2,10,'k','filled')
            axis square
            
            xlim([0 1])
            ylim(xlim)
            
            [r p]=corr(v1,v2,'rows','complete');
            text(0.8,0.1,num2str(r))
            
            m=m+1;
            
            
        end
        
    end
end




end


