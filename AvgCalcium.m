%% This function/script get the raw time series data from the results object. It removes the nan values and translate the intensities values into calcium concentraiton by interpolating with the baseline flou4 intensity. 
%{
  Inputs:
         Tca: vector of time points
         
         Ca: matrix of intensity time series. It is assumed that the rows
         are the timepoints and columns are the cell numbers
         
         baseInd: a vector of indices that indicate the range of baseline intensity values. 
         outInd: the vector of indices that indicate the data points that should be discarded.
         baseInten: the corresponding species concentration for the
         baseline intensity(in units of micromolar)
  Outputs:
%}
%% Get calcium data
[Ca,Tca] = R.getTimeseriesData('Ca','B10');
%%
figure(2);clf; prm = randperm(size(Ca,2)); plot(Tca,Ca(:,prm(1:10)));
title('Sample individual time series of B10 (0uM)');
xlabel('Time (sec)');
ylabel('Raw Flou4 intensity');
%%
baseInd = 1:10;
outInd = 11:17;
baseInten = 50;

%%
%function [TcaAvg,CaAvg] = AvgCalcium(Tca,Ca,baseInd,outInd, baseInten)    
    CaTemp = Ca;
    %remove the cells with nan entries
    [row,col] = find(isnan(CaTemp)); CaTemp(:,unique(col)) = [];
    % calculate the vectors of baseline intensity level   
    baseLevel = mean(CaTemp(baseInd,:),1);     
    % remove the outlier timepoints
    TcaTemp = Tca;
    CaTemp(outInd,:) =nan;
    %TcaTemp(outInd) =[];

    % divide the whole Ca matrix by the baseLevel
    for i = 1:size(CaTemp,1)
        CaTemp(i,:) = CaTemp(i,:)./baseLevel;
    end
    % Multiply the whole Ca matrix by the baseInten to convert the
    % proportional increase of intensity to species concentration
    CaTemp = CaTemp.*baseInten;
    % Get the average behavior by averaging the time series of all the
    % cells
    CaAvg0 = mean(CaTemp,2);
    TcaAvg0 = TcaTemp;
    %%

    %CaArray = {CaAvg10; CaAvg3_3; CaAvg1; CaAvg033; CaAvg01; CaAvg0};
    %TcaArray ={TcaAvg10; TcaAvg3_3; TcaAvg1; TcaAvg033; TcaAvg01; TcaAvg0};
    CaArray = {CaAvg10; CaAvg1; CaAvg01; CaAvg0};
    TcaArray ={TcaAvg10; TcaAvg1; TcaAvg01; TcaAvg0};
    %%
    % new figure    
    figure(3);
    clr = jet(length(CaArray));
    for i = 1:length(CaArray)
        plot(TcaArray{i},CaArray{i},'color',clr(i,:),'LineWidth',2);
        hold on;
    end
    %legend('10uM','3.33uM','1uM','0.33uM','0.1uM','Control');
    legend('10uM','1uM','0.1uM','Control');
    h2 = get(gcf,'CurrentAxes');
    title('Average Calcium Response 12/14/2013 run 2');
    h2 = get(gca,'Title');
    set(h2,'FontSize', 20);
    xlabel('Time sec');    
    h2 = get(gca,'XLabel');
    set(h2,'FontSize',20);
    ylabel('Concentration (micromolar)');    
    h2 = get(gca,'YLabel');
    set(h2,'FontSize',20);
    %set(gca,'LineWidth',3);