%% Init
clear
close all
clc
stkshow('close','all');

% place where the copy of the EKARevAnalysisFunctions are (the EKARev repo)
addpath /home/rwollman/Projects/InVitroWounding/EKARevAnalysisFunctions/

%% Key user input - the path where the images are
pth='/data/Images/anna/InVitroWounding/EKARev_ATP_GM6001_2013Jul09'; 


%% Init MD / R and create the Props per position
MD=Metadata(pth);
R = MultiPositionSingleCellResults(pth); 
R.PosNames = unique(MD,'Position');  

% Props are the property that 
Props =  MD.NewTypes; 

%% Main loop
t0=now; 
for i=1:numel(R.PosNames)
    % output timing
    fprintf('%s %s\n',R.PosNames{i},datestr(now-t0,13)); 
    
    % trick to prevent segmenting if segmentation exist already
    Lbl = @() segmenetEKARev(MD,R.PosNames{i}); 
    Lbl = setLbl(R,Lbl,R.PosNames{i}); 
    
    % measure intensities
    [Erk,Terk,FretStk] = measureEKARevFRET(MD,R.PosNames{i},Lbl);
    
    % add measurements
    addTimeSeries(R,'Erk',Erk,Terk,R.PosNames{i}); 
    add(R,'FretStack',FretStk,R.PosNames{i}); 
    
    % add Position properties
    for j=1:numel(Props)
        tmp=unique(MD,Props{j},'Position',R.PosNames{i});
        setProperty(R,Props{j},tmp(1),R.PosNames{i});
    end
        
end

%% save results
R.analysisScript=fullfile(pwd,mfilename); 
R.saveResults; 
