function []=dot_and_line_plot(inputArray)

toPlot1=nan(length(inputArray),1);
toPlot2=toPlot1;

for i=1:length(toPlot1)
    toPlot1(i)=mean(inputArray{i},'omitnan');
    toPlot2(i)=std(inputArray{i},[],'omitnan')./sqrt(length(inputArray{i}));
end

hold on
errorbar(1:length(toPlot1),toPlot1,toPlot2,'.k','CapSize',0)
scatter(1:length(toPlot1),toPlot1,50,'k','filled')
xlim([0.5 length(toPlot1)+0.5])
ylim([0 Inf])

end