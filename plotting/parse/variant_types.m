%function that returns variant type distribution
function [distOut,arrayOut]=variantTypes4(inputArray)

%protein-altering
variantArray{1}={'disruptive_inframe_deletion','disruptive_inframe_insertion',...
    'frameshift_variant','frameshift_variant&start_lost',...
    'frameshift_variant&stop_gained','frameshift_variant&stop_lost&splice_region_variant',...
    'inframe_deletion','inframe_insertion','initiator_codon_variant',...
    'intergenic_region','missense_variant','start_lost',...
    'start_lost&inframe_insertion','stop_gained','stop_gained&disruptive_inframe_deletion',...
    'stop_gained&inframe_insertion','stop_lost&splice_region_variant'};
%syn
variantArray{2}={'synonymous_variant'};

%intron and splice
% variantArray{3}={'intron_variant','non_coding_transcript_variant',...
%     'splice_donor_variant&intron_variant','splice_region_variant&intron_variant',...
%     'splice_region_variant&stop_retained_variant'};
%UTR --  has reads in mRNAseq
variantArray{3}={'intron_variant','non_coding_transcript_variant',...
    'splice_donor_variant&intron_variant','splice_region_variant&intron_variant',...
    'splice_region_variant&stop_retained_variant','UTR','intragenic_variant','intergenic_region'};
%intergenic
%variantArray{4}={'intergenic_region'};

arrayOut=zeros(1,length(inputArray));

for i=1:length(variantArray)

    arrayOut(ismember(inputArray,variantArray{i}))=i;

end

%arrayOut(arrayOut==0)=length(variantArray)+1;

%for i=1:length(variantArray)+1
for i=1:length(variantArray)
     
    distOut(i)=sum(arrayOut==i);
    
end

distOut=distOut./sum(distOut);

