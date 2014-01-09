%{
This script allows for the user to quickly inspect the raw data from the individual wells 
%}
%% define the directory of the raw data, its metadata, and the well to be analyzed
pth = '/data2/Images/Jason/CalciumDynamics/ATPdoseResponse_2013Dec18';
md = Metadata(pth);
pos = {'B04'};
%% read the nuclei stack
nuc = stkread(md,'Position',pos,'Channel','DeepBlue');
%% read the calcium stack
ca = stkread(md,'Position',pos,'Channel','Yellow');
%% call the function to process the data
R = ProcessCalciumData(pth,pos);
%% print out sample invidivual time series
[Ca,Tca] = R.getTimeseriesData('Ca',pos{:});
figure; clf; prm = randperm(size(Ca,2)); plot(Tca,Ca(:,prm(1:10)));
title('Sample individual time series');
xlabel('Time (sec)');
ylabel('Raw Flou4 intensity');
legend(num2str(prm(1)),num2str(prm(2)),num2str(prm(3)),num2str(prm(4)),num2str(prm(5)),num2str(prm(6)),num2str(prm(7)), num2str(prm(8)),num2str(prm(9)),num2str(prm(10)));