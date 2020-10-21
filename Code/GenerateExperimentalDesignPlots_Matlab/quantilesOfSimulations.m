%% Plotting the increasing gamma and increasing d designs
% Date : 2020-10-20
% Author: Johannes Borgqvist
% Description:
% We generate tex-files which can subsequently
% be plotted in pgfplot.You need to give the 
% folders from which you read the data and
% provide the destination folders where the
% plotting is made. Then, the program will 
% make the tex-files which is the quantiles of 
% the two designs increasing gamma and increasing d. 
clear,clc

%% INCREASING GAMMA

% We start by reading the variable vector
data = readtable('../../Results/increasingGamma/Classical/1/gammaVec.csv');
temp = table2array(data);
gammaVec = temp(2:end,2);
extraOption = false; % whether to forget plot or not in the legend
% Define the number of experiments
%nuOfExperiments = length(gammaVec);
%% -----------------------------------------------------------------------------------
% NUMBER OF POLES: CLASSIC
% Classic
str = '../../Results/increasingGamma/Classical/nuOfPoles_increasingGamma_Classic.xlsx';
nuOfPoles = xlsread(str);
strDest = '../../Figures/increasing_gamma/nuOfPoles/Input/';
destinationFolder = {[strDest, 'nuOfPoles_Classical.tex'], [strDest, 'nuOfPoles_Classical_High.tex'], [strDest, 'nuOfPoles_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'nuOfPolesClassical';
generatePlots(nuOfPoles,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);
% Non-Classic
str = '../../Results/increasingGamma/NonClassical/nuOfPoles_increasingGamma_NonClassic.xlsx';
nuOfPoles = xlsread(str);
strDest = '../../Figures/increasing_gamma/nuOfPoles/Input/';
destinationFolder = {[strDest, 'nuOfPoles_NonClassical.tex'], [strDest, 'nuOfPoles_NonClassical_High.tex'], [strDest, 'nuOfPoles_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'nuOfPolesNonClassical';
generatePlots(nuOfPoles,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);
%% -----------------------------------------------------------------------------------
% UMAX AND UMIN
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CLASSIC UMIN AND UMAX
% The number of data points are
nuOfRepititions = 20;
% Read data: umin
dataStr = '/uMin.csv';
Case = 'increasingGamma/Classical/';
uMinClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: umin
strDest = '../../Figures/increasing_gamma/maxMin/Input/';
destinationFolder = {[strDest, 'uMin_Classical.tex'], [strDest, 'uMin_Classical_High.tex'], [strDest, 'uMin_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'uMinClassical';
generatePlots(uMinClassic,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);
% umax we forget in the legend
extraOption = true;
% Read data: umax
dataStr = '/uMax.csv';
Case = 'increasingGamma/Classical/';
uMaxClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: umax
strDest = '../../Figures/increasing_gamma/maxMin/Input/';
destinationFolder = {[strDest, 'uMax_Classical.tex'], [strDest, 'uMax_Classical_High.tex'], [strDest, 'uMax_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'uMaxClassical';
generatePlots(uMaxClassic,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NON-CLASSIC UMIN AND MAX
% umax we do not forget in the legend
extraOption = false;
% The number of data points are
nuOfRepititions = 19;
% Read data: umin
dataStr = '/uMin.csv';
Case = 'increasingGamma/NonClassical/';
uMinNonClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: umin
strDest = '../../Figures/increasing_gamma/maxMin/Input/';
destinationFolder = {[strDest, 'uMin_NonClassical.tex'], [strDest, 'uMin_NonClassical_High.tex'], [strDest, 'uMin_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'uMinNonClassical';
generatePlots(uMinNonClassic,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);
% umax we forget in the legend
extraOption = true;
% Read data: umax
dataStr = '/uMax.csv';
Case = 'increasingGamma/NonClassical/';
uMaxNonClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: umax
strDest = '../../Figures/increasing_gamma/maxMin/Input/';
destinationFolder = {[strDest, 'uMax_NonClassical.tex'], [strDest, 'uMax_NonClassical_High.tex'], [strDest, 'uMax_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'uMaxNonClassical';
generatePlots(uMaxNonClassic,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);
% for tPol we forget again
extraOption = false;

%% -----------------------------------------------------------------------------------
% POLARISATION TIME
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CLASSIC POLARISATION TIME
% The number of data points are
nuOfRepititions = 20;
% Read data:  
dataStr = '/tPole.csv';
Case = 'increasingGamma/Classical/';
tPoleClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: tPole
strDest = '../../Figures/increasing_gamma/polTime/Input/';
destinationFolder = {[strDest, 'polTime_Classical.tex'], [strDest, 'polTime_Classical_High.tex'], [strDest, 'polTime_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'polTimeClassical';
generatePlots(tPoleClassic,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NON-CLASSIC POLARISATION TIME
% The number of data points are
nuOfRepititions = 19;
%Read data:  
dataStr = '/tPole.csv';
Case = 'increasingGamma/NonClassical/';
tPoleNonClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
%Write the data: tPole
strDest = '../../Figures/increasing_gamma/polTime/Input/';
destinationFolder = {[strDest, 'polTime_NonClassical.tex'], [strDest, 'polTime_NonClassical_High.tex'], [strDest, 'polTime_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'polTimeNonClassical';
generatePlots(tPoleNonClassic,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);


%% -----------------------------------------------------------------------------------
% POLE RATIO
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CLASSIC POLE RATIO
% The number of data points are
nuOfRepititions = 20;
% Read data:  
dataStr = '/ratioPole.csv';
Case = 'increasingGamma/Classical/';
ratioPoleClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: ratioPole
strDest = '../../Figures/increasing_gamma/ratioPole/Input/';
destinationFolder = {[strDest, 'ratioPole_Classical.tex'], [strDest, 'ratioPole_Classical_High.tex'], [strDest, 'ratioPole_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'ratioPoleClassical';
generatePlots(ratioPoleClassic,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% NON-CLASSIC POLE RATIO
% The number of data points are
nuOfRepititions = 19;
% Read data:  
dataStr = '/ratioPole.csv';
Case = 'increasingGamma/NonClassical/';
ratioPoleNonClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: ratioPole
strDest = '../../Figures/increasing_gamma/ratioPole/Input/';
destinationFolder = {[strDest, 'ratioPole_NonClassical.tex'], [strDest, 'ratioPole_NonClassical_High.tex'], [strDest, 'ratioPole_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'ratioPoleNonClassical';
generatePlots(ratioPoleNonClassic,destinationFolder,gammaVec,colourStr,fillSpaceStr,opacity,extraOption);



%% INCREASING d

% We start by reading the variable vector
data = readtable('../../Results/increasing_d/Classical/1/dVec.csv');
temp = table2array(data);
dVec = temp(2:end,2);
% Define the number of experiments
nuOfExperiments = length(dVec);
% The number of data points are
nuOfRepititions = 20;
%% -----------------------------------------------------------------------------------
% UMAX AND UMIN
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CLASSIC UMIN AND UMAX
% Read data: umin
dataStr = '/uMin.csv';
Case = 'increasing_d/Classical/';
uMinClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: umin
strDest = '../../Figures/increasing_d/maxMin/Input/';
destinationFolder = {[strDest, 'uMin_Classical.tex'], [strDest, 'uMin_Classical_High.tex'], [strDest, 'uMin_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'uMinClassical';
generatePlots(uMinClassic,destinationFolder,dVec,colourStr,fillSpaceStr,opacity,extraOption);
% No legend for uMax
extraOption = true;
% Read data: umax
dataStr = '/uMax.csv';
Case = 'increasing_d/Classical/';
uMaxClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: umax
strDest = '../../Figures/increasing_d/maxMin/Input/';
destinationFolder = {[strDest, 'uMax_Classical.tex'], [strDest, 'uMax_Classical_High.tex'], [strDest, 'uMax_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'uMaxClassical';
generatePlots(uMaxClassic,destinationFolder,dVec,colourStr,fillSpaceStr,opacity,extraOption);
 % We remember uMin
extraOption = false;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%NON-CLASSIC UMIN AND MAX
%CLASSIC UMIN AND UMAX
% Read data: umin
dataStr = '/uMin.csv';
Case = 'increasing_d/NonClassical/';
uMinNonClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: umin
strDest = '../../Figures/increasing_d/maxMin/Input/';
destinationFolder = {[strDest, 'uMin_NonClassical.tex'], [strDest, 'uMin_NonClassical_High.tex'], [strDest, 'uMin_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'uMinNonClassical';
generatePlots(uMinNonClassic,destinationFolder,dVec,colourStr,fillSpaceStr,opacity,extraOption);
% uMax we forget
extraOption = true;
% Read data: umax
dataStr = '/uMax.csv';
Case = 'increasing_d/NonClassical/';
uMaxNonClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: umax
strDest = '../../Figures/increasing_d/maxMin/Input/';
destinationFolder = {[strDest, 'uMax_NonClassical.tex'], [strDest, 'uMax_NonClassical_High.tex'], [strDest, 'uMax_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'uMaxNonClassical';
generatePlots(uMaxNonClassic,destinationFolder,dVec,colourStr,fillSpaceStr,opacity,extraOption);
% tPol we have in the legend
extraOption = false;
%% -----------------------------------------------------------------------------------
% POLARISATION TIME
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CLASSIC POLARISATION TIME
% Read data:  
dataStr = '/tPole.csv';
Case = 'increasing_d/Classical/';
tPoleClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: tPole
strDest = '../../Figures/increasing_d/polTime/Input/';
destinationFolder = {[strDest, 'polTime_Classical.tex'], [strDest, 'polTime_Classical_High.tex'], [strDest, 'polTime_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'polTimeClassical';
generatePlots(tPoleClassic,destinationFolder,dVec,colourStr,fillSpaceStr,opacity,extraOption);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%NON-CLASSIC POLARISATION TIME
%Read data:  
dataStr = '/tPole.csv';
Case = 'increasing_d/NonClassical/';
tPoleNonClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
%Write the data: tPole
strDest = '../../Figures/increasing_d/polTime/Input/';
destinationFolder = {[strDest, 'polTime_NonClassical.tex'], [strDest, 'polTime_NonClassical_High.tex'], [strDest, 'polTime_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'polTimeNonClassical';
generatePlots(tPoleNonClassic,destinationFolder,dVec,colourStr,fillSpaceStr,opacity,extraOption);


% -----------------------------------------------------------------------------------
%POLE RATIO
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CLASSIC POLE RATIO
% Read data:  
dataStr = '/ratioPole.csv';
Case = 'increasing_d/Classical/';
ratioPoleClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: ratioPole
strDest = '../../Figures/increasing_d/ratioPole/Input/';
destinationFolder = {[strDest, 'ratioPole_Classical.tex'], [strDest, 'ratioPole_Classical_High.tex'], [strDest, 'ratioPole_Classical_Low.tex']};
colourStr = 'blue'; 
opacity = 0.1;
fillSpaceStr = 'ratioPoleClassical';
generatePlots(ratioPoleClassic,destinationFolder,dVec,colourStr,fillSpaceStr,opacity,extraOption);

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%NON-CLASSIC POLE RATIO
%Read data:  
dataStr = '/ratioPole.csv';
Case = 'increasing_d/NonClassical/';
ratioPoleNonClassic = readData(dataStr,Case,nuOfRepititions); % Read the data
% Write the data: ratioPole
strDest = '../../Figures/increasing_d/ratioPole/Input/';
destinationFolder = {[strDest, 'ratioPole_NonClassical.tex'], [strDest, 'ratioPole_NonClassical_High.tex'], [strDest, 'ratioPole_NonClassical_Low.tex']};
colourStr = 'orange'; 
opacity = 0.1;
fillSpaceStr = 'ratioPoleNonClassical';
generatePlots(ratioPoleNonClassic,destinationFolder,dVec,colourStr,fillSpaceStr,opacity,extraOption);
