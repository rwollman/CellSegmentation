%{
This script allows for the user to quickly inspect the raw data from the individual wells 
%}
%% Init

clear
close all
clc
stkshow('close','all')

%%
%posArray = {'B04','B05','B06','B08','B09','B10'};
posArray = {'B03'};
TcaArray = cell(length(posArray),1);
meanCaArray = cell(length(posArray),1);
for arrayInd = 1:length(posArray);
    disp(strcat('Analyzing well ',posArray{arrayInd}));
    %% define the directory of the raw data, its metadata, and the well to be analyzed
    pth = '/data2/Images/Jason/CalciumDynamics/ATPdoseResponse_2014Mar03';
    md = Metadata(pth);
    pos = posArray{arrayInd};

    %% read the nuclei stack
    %{
    nuc = stkread(md,'Position',pos,'Channel','DeepBlue');

    %% read the calcium stack
    ca = stkread(md,'Position',pos,'Channel','Yellow','resize',0.25);
    %}

    %%
    [Lbl,NucLabels,T,msk] = segmentNucleiOnly(md,pos);

    %%
    [Ca,T,CaStk] = measureStackIntensities(md,pos,Lbl,'channel','Yellow','background',false,'register',true); 
    %% get the Time vector
    Tca = md.getSpecificMetadata('TimestampFrame','Channel','Yellow','Position',pos);
    Tca = cat(1,Tca{:});
    [Tca ,ordr]= sort(Tca);
    Tca = (Tca - Tca(1))*86400;
    TcaArray{arrayInd} = Tca;
    % get the normalized calcium time series
    meanCa = mean(Ca,2);
    meanCa  = (meanCa/(mean(meanCa(1:7))))*50;
    meanCaArray{arrayInd} = meanCa;
    
end
%%
%{
clr = jet(length(posArray));
figure(1);
clf;
hold on;
for i =1:size(clr,1)
    plot( TcaArray{i},meanCaArray{i},'color',clr(i,:));
end
legend('10 uM','3.33 uM','1 uM', '0.33 uM', '0.1 uM', ' control');
title('Experimental Data');
xlabel('Time (seconds)');
ylabel('concentration (micromolar)');
%}
%{
%% plot all the time series
TcaArray = { TcaB04, TcaB05, TcaB06, TcaB08, TcaB09, TcaB10};
meanCaArray ={meanCaB04, meanCaB05,meanCaB06, meanCaB08, meanCaB09, meanCaB10 }; 
clr = jet(6);
figure(1);
clf;
hold on;
for i = 1:size(clr,1)
    plot(TcaArray{i},meanCaArray{i},'color',clr(i,:));
end
legend('10 uM','3.33 uM','1 uM', '0.33 uM', '0.1 uM', ' control');

%}
%%
%{
R = MultiPositionSingleCellResults(pth); % ,'timefunc',@(t) t<datenum('18-Dec-2013 12:20:30')
setLbl(R,Lbl,R.PosNames{j});
addTimeSeries(R,'Ca',Ca,Tca,R.PosNames{j});
%}
% R = ProcessCalciumData(pth,pos,'timefunc',@(t) t<datenum('18-Dec-2013 12:20:30'));
% 
% % print out sample invidivual time series
% [Ca,Tca] = R.getTimeseriesData('Ca',pos{:});
% figure; clf; prm = randperm(size(Ca,2)); plot(Tca,Ca(:,prm(1:10)));
% title('Sample individual time series');
% xlabel('Time (sec)');
% ylabel('Raw Flou4 intensity');
% legend(num2str(prm(1)),num2str(prm(2)),num2str(prm(3)),num2str(prm(4)),num2str(prm(5)),num2str(prm(6)),num2str(prm(7)), num2str(prm(8)),num2str(prm(9)),num2str(prm(10)));