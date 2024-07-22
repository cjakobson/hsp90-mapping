function v_dist = calculate_telomere_distance(dependency_directory,output_directory)


variant_input=readtable([dependency_directory 'variantInfoStructure.csv']);

telo_input=readtable([dependency_directory 'chromosome_length.txt']);


v_dist=nan(height(variant_input),1);

for i=1:length(v_dist)
   
    temp_chr=variant_input.chr(i);
    temp_pos=variant_input.pos(i);
    
    temp_chr_length=telo_input.Var3(telo_input.Var1==temp_chr);
    
    if ~isempty(temp_chr_length)
        temp_dist=abs(temp_chr_length-temp_pos);

        v_dist(i)=min([temp_pos,temp_dist]);
        %normalize to chromosome length
        v_dist(i)=v_dist(i)/temp_chr_length;
    end
    
end


end


