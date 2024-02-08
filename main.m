% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Tianwei
% Date: 06/01/2024
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add to the path
fullFileName = mfilename('fullpath');
addpath(genpath(fileparts(fullFileName)));


%% Fig 1
Fig1bc = plotHandtrj_emg;


%% Fig 2
Fig2 = Exam_Cell; %pick neuron data in the exampleN folder, plot three psths aligned to GO MO ME

%% Fig 3
Fig3 = plotFig3; %populaiton fit performance on in 3 coordinations; an example fit result of neuron 142

%% Fig 4
Fig4 = plotTuningCurve;

%% Fig 5
Fig5a = plot_perm_Carray; %perm_dir
Fig5b = plot_perm_G; %data_Guass

%% Fig 6
Fig6ac = PCALDA_NC;
Fig6b = deco;

%% Fig 7
Fig7 = Simulate_PV;






