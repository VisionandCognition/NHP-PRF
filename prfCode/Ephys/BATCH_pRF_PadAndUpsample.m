warning off;clear all; clc;

addpath(genpath('../../prfCode'));

fprintf('=================================================\n');
fprintf('Padding & upsampling for Danny \n');
fprintf('=================================================\n');

pRF_prepdata_avg('pRF_PrepDatalist_Danny',1)

fprintf('\n=================================================\n');
fprintf('Padding & upsampling for Eddy \n');
fprintf('=================================================\n')

pRF_prepdata_avg('pRF_PrepDatalist_Eddy',1)

fprintf('== ALL DONE! ==\n')

rmpath(genpath('../../prfCode'));