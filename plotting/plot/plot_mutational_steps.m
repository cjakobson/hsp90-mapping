function [] = plot_mutational_steps(loci_to_plot,ref_allele,condition_to_plot,...
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

%calculate pVals between genotypes
p_mat=nan(length(to_plot));
for i=1:length(to_plot)
    
    v_fit(i)=mean(to_plot{i},'omitnan');
    
    for j=1:length(to_plot)
        
        [h,p_mat(i,j)]=ttest2(to_plot{i},to_plot{j});
        
    end
    
end

allowed_mat=p_mat>0.05;

%also fitness for line widths
v_fit=v_fit+1.5;



  


%make mutational step plots
for i=1:length(to_plot)
    v_mean(i)=mean(to_plot{i},'omitnan');
end

%connect relevant dots
temp_mean=v_mean(1:2:end);
temp_dist=v_dist(1:2:end);

m=1;
for i=0:3
    
    idx1=find(temp_dist==i);
    
    for j=i+1
        
        idx2=find(temp_dist==j);
        
        for k=1:length(idx1)
            
            for l=1:length(idx2)
                
                v_delta{1}(m)=temp_mean(idx2(l))-temp_mean(idx1(k));
                m=m+1;
                
            end
            
        end
        
    end
end


temp_mean=v_mean(2:2:end);
temp_dist=v_dist(2:2:end);

m=1;
for i=0:3
    
    idx1=find(temp_dist==i);
    
    for j=i+1
        
        idx2=find(temp_dist==j);
        
        for k=1:length(idx1)
            
            for l=1:length(idx2)
                
                v_delta{2}(m)=temp_mean(idx2(l))-temp_mean(idx1(k));
                m=m+1;
                
            end
            
        end
        
    end
end


hold on
easy_box_with_dots(v_delta)
ylim([-2 1])
title(condition_to_plot)
[h p]=ttest2(to_plot{1},to_plot{1});
text(1.5,0.5,num2str(p))



end


