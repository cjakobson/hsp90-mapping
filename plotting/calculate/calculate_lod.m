function [lod_1D,perm_thresh] = calculate_lod(trait_to_use,n_perms,dependency_directory,output_directory)


model_genotypes = parse_genotypes(dependency_directory,output_directory);


%coarse LOD scoring to start (no adjustments) to ID hotspots
lod_1D=nan(1,length(model_genotypes(1,:)));
for q = 1:length(model_genotypes(1,:))
    r = corr(model_genotypes(:,q),trait_to_use,'rows','complete');
    lod_1D(q) = -length(trait_to_use)*log(1-r^2)/(2*log(10));
end

%permute N times
perm_lod=[];
for i=1:n_perms
    
    temp_lod=nan(1,length(model_genotypes(1,:)));
    
    rng(i)
    
    temp_trait=trait_to_use(randperm(length(trait_to_use)));
    
    for q = 1:length(model_genotypes(1,:))
        r = corr(model_genotypes(:,q),temp_trait,'rows','complete');
        temp_lod(q) = -length(temp_trait)*log(1-r^2)/(2*log(10));
    end
    
    perm_lod=[perm_lod temp_lod];
    
end

sorted_perm=sort(perm_lod,'descend');

%use 5% FDR
desired_fdr=0.05;
perm_thresh=sorted_perm(ceil(length(sorted_perm)*desired_fdr*n_perms));



end


