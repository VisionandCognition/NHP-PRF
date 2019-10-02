function pRF_prepCrossVal(monkey)

clc;
startfld=pwd;
fprintf(['Prepping crossvalidation data input for ' monkey '\n']);

cd ..; cd Data; cd MRI; cd avg; cd(monkey);

R = load('AllSessions-avg-odd.mat');
R(2) = load('AllSessions-avg-even.mat');

sess_meanBOLD = {R(1).sess_meanBOLD,R(2).sess_meanBOLD};
sess_wmeanBOLD = {R(1).sess_wmeanBOLD,R(2).sess_wmeanBOLD};
sess_medianBOLD = {R(1).sess_medianBOLD,R(2).sess_medianBOLD};
sess_sdBOLD = {R(1).sess_sdBOLD,R(2).sess_sdBOLD};

sess_meanBOLD_inv = {R(1).sess_meanBOLD_inv,R(2).sess_meanBOLD_inv};
sess_wmeanBOLD_inv = {R(1).sess_wmeanBOLD_inv,R(2).sess_wmeanBOLD_inv};
sess_medianBOLD_inv = {R(1).sess_medianBOLD_inv,R(2).sess_medianBOLD_inv};
sess_sdBOLD_inv = {R(1).sess_sdBOLD_inv,R(2).sess_sdBOLD_inv};

stim.norm = {R(1).stim.norm, R(2).stim.norm};
stim.inv = {R(1).stim.inv, R(2).stim.inv};

clear R;

cd(startfld)
[~,~,~] = mkdir(fullfile('..','.Data','MRI','cv',MONKEY{m}));
save('AllSessions-avg-cv.mat');
cd(startfld)
