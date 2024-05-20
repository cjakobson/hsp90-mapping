function []=plot_allele_combintations(loci_to_plot,ref_allele,condition_to_plot,...
    dependency_directory,output_directory)


set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



load([dependency_directory 'radFilename.mat'])
load([dependency_directory 'radTrait.mat'])

trait_idx1=ismember(filename,[condition_to_plot '-rad']);
trait_idx2=ismember(filename,[condition_to_plot '+rad']);

model_genotypes = parse_genotypes(dependency_directory,output_directory);

for i=1:length(loci_to_plot)
    
    temp_locus=loci_to_plot(i);
    
    v_genotype{i,1}=model_genotypes(:,temp_locus)==ref_allele(i);
    v_genotype{i,2}=model_genotypes(:,temp_locus)==-ref_allele(i);
    
end


m=1;
for i=1:2
    for j=1:2
        for k=1:2
            %-/+rad
            v_temp=logical(v_genotype{1,i}.*v_genotype{2,j}.*v_genotype{3,k});
            to_plot{m}=trait{trait_idx1}(v_temp);
            
            %mutational distance
            v_dist(m)=sum([i==2 j==2 k==2]);
            
            m=m+1;
            
            
            to_plot{m}=trait{trait_idx2}(v_temp);
            
            %mutational distance
            v_dist(m)=sum([i==2 j==2 k==2]);
            
            m=m+1;
            
        end
    end
end



hold on
easyBox(to_plot)
ylim([-2 2])
title(condition_to_plot)




end



