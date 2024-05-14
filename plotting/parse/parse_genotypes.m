function model_genotypes = parse_genotypes(dependency_directory,output_directory)

    
    
load([dependency_directory 'phasedVLCgenotype.mat'])
genotypes=phasedVLCgenotype;
clear phasedVLCgenotype


load([dependency_directory 'radTrait.mat'])

n_strains=length(trait{1});

%truncate to measured strains
genotypes=genotypes(1:n_strains,:);

%MODEL A
%ignore hets; homoRM gets 1; homoYJM gets -1
%final: only nStrains columns
[~,temp]=size(genotypes);
n_loci=temp/4;    

model_genotypes=zeros(n_strains,n_loci);
for i=1:n_strains
    v_genotype=genotypes(i,1:n_loci)-genotypes(i,(n_loci+1):(2*n_loci));
    model_genotypes(i,:)=v_genotype;
end

clear genotypes



end
