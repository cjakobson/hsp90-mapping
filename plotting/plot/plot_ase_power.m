function [] = plot_ase_power(dependency_directory,output_directory)

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



%calculate what AF yields q<0.05 as f(read depth)
%median true depth is 60 reads/nt
%assume 8000 tests for bonferroni
v_maf=0.5:0.01:1;
v_reads=floor(10.^(1:0.1:4));
for i=1:length(v_maf)%underlying MAF
    
    for j=1:length(v_reads)
        
        q_mat(i,j)=(1-binocdf(v_maf(i)*v_reads(j),v_reads(j),0.5)+1e-17)*8000;
        
    end
    
end



hold on
plot(-log10(q_mat(:,v_reads==100)))
plot(-log10(q_mat(:,v_reads==1000)))
ylabel('-log_{10}qVal')
xlabel('true MAF')
xticks(1:10:length(v_maf))
xticklabels(v_maf(1:10:length(v_maf)))
xtickangle(45)
axis square

legend({'100 reads/nt','1000 reads/nt'},'Location','Southeast')

end


