%{
This script allows for the user to quickly inspect the raw data from the individual wells 
%}
%% define the directory of the raw data, its metadata, and the well to be analyzed
pth = '/data2/Images/Jason/CalciumDynamics/ATPdoseResponse_2013Dec18';
md = Metadata(pth);
pos = 'B04';
%% read the nuclei stack
nuc = stkread(md,'Position',pos,'Channel','DeepBlue');
%% read the calcium stack
ca = stkread(md,'Position',pos,'Channel','Yellow');
