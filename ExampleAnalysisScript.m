%% Init
clear
close all
clc
stkshow('close','all');

%% Key user input - the path where the images are
pth='/data/Images/where/your/dataset/is'; 


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
    [Ca,Tca,CaStk] = measureMCF10ACa2(MD,grp{i},Lbl); 
    
    % create movies
    R.stack2movie(FretStk,sprintf('%s_Erk',R.PosNames{i}));
    R.stack2movie(CaStk,sprintf('%s_Erk',R.PosNames{i}),'colormap',gray(256)); 
    
    % add measurements
    addTimeSeries(R,'Erk',Erk,Terk,R.PosNames{i}); 
    addTimeSeries(R,'Ca',Ca,Tca,R.PosNames{i}); 
    
    % add Position properties
    for j=1:numel(Props)
        tmp=unique(MD,Props{j},'Position',R.PosNames{i});
        setProperty(R,Props{j},tmp(1),R.PosNames{i});
    end
        
end

%% save results
R.analysisScript=fullfile(pwd,mfilename); 
R.saveResults; 
