%prepare input data for Hsp90 mapping

clear

tic

set(0,'DefaultLineLineWidth',1)
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultAxesLineWidth',1)

blue=[43 172 226]./256;
orange=[248 149 33]./256;
grey=[128 128 128]./256;



filebase='/Users/christopherjakobson/';


code_directory=[filebase 'Documents/GitHub/hsp90-mapping/'];
dependency_directory=[filebase '/Dropbox/JaroszLab/hsp90mapping/hsp90-mapping-dependencies/'];
output_directory=[filebase 'Dropbox/JaroszLab/hsp90mapping/manuscript-plots/'];

addpath([code_directory 'plotting'])
addpath([code_directory 'plotting/parse'])
addpath([code_directory 'plotting/calculate'])
addpath([code_directory 'plotting/plot'])
addpath([code_directory 'data-prep'])



analyze_phenotyping(dependency_directory,output_directory)




toc




