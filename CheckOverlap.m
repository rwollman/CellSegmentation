%%This script get the overlapping labels that do not have nans
[Ca,Tca] = R.getTimeseriesData('Ca','B07');
Lbl = R.getLbl('B07');
%% show all the cells with nans
[row,col] = find(isnan(Ca));
troubleMaker = unique(col);
label = Lbl.Lbl;
nanLabel = max(max(max(label))) +1;
% change all the cells with nans to a new cell id which is 100 more than
% the existing maximum cell id
for i =1:size(label,3)
    nanPresent = ismember(label(:,:,i), troubleMaker);
    labeli = label(:,:,i);
    labeli(nanPresent) = nanLabel;
    label(:,:,i) = labeli;
end
figure(1); imshow(label(:,:,1) == nanLabel);
%% Get the cell id's that overlap with each other
bwArray = cell(size(Lbl.Lbl,3),1);
for i =1:length(bwArray)
    bwArray{i} = logical(Lbl.Lbl(:,:,i));
end
bwArray = cat(3,bwArray{:});

%% Get the overlap of labels between each of the images
overlap = cell(size(Lbl.Lbl,3)-1,1);
for i = 1:length(overlap)
    overlap{i} = bwArray(:,:,1) & bwArray(:,:,i+1);
end
overlap = cat(3,overlap{:});
%% Get the overlap pairs
pairs = cell(size(Lbl.Lbl,3)-1,1);
lbl1 = Lbl.Lbl(:,:,1);
for i=1:length(pairs)
    lblComp = Lbl.Lbl(:,:,i+1); 
    pairs{i} =[lbl1(overlap(:,:,i)), lblComp(overlap(:,:,i))];
    pairs{i} = unique(pairs{i},'rows');
    p = pairs{i};
end
%% Get the overlap pairs into sparse matrix
S = cell(length(pairs),1);
for i = 1:length(S)
    p = pairs{i};
    S{i} = sparse(p(:,1),p(:,2),1);
end
%% Get the volumes of cells corresponding to cell id's from the labeled images
%cellVolume is a matrix that keeps track of the volumes of cells
labelSet = unique(label(:,:,1));
cellVolumeMat = zeros(length(labelSet),size(label,3)+1);
cellVolumeMat(:,1) = labelSet;
for i = 1:size(label,3)
    currLabel = label(:,:,i);
    countCell = hist(currLabel(:),labelSet);
    cellVolumeMat(:,i+1) = countCell';
end

cellVolume = cell(size(label,3),1);
for i = 1:size(cellVolume,1)
    currLabel = label(:,:,i);
    cellSet = unique(currLabel(:));
    countCell = hist(currLabel(:),cellSet);
    cellVolume{i} = [cellSet, countCell'];
end

%% Obtain information on how much cells change to other cells from one labeled image to another
% cellOverlap keeps track of the change in cell volumes from one label to
% another
cellOverlap = cell(size(label,3)-1,1);
for i =1:size(cellOverlap,1)
    oldLabel = label(:,:,i);
    currLabel  = label(:,:,i+1);
    matrixOverlap = bwArray(:,:,i) & bwArray(:,:,i+1);
    currLabelPair = [oldLabel(matrixOverlap) currLabel(matrixOverlap)];
    % remove entries of zero from currLabelPair
    [row ,col] = find(currLabelPair==0);
    currLabelPair(row,:) = [];
    %labelPair{i} = [oldLabel(matrixOverlap) currLabel(matrixOverlap)]; 
    % count the number of each of the overlap occurrences 
    [C,ia,ic] = unique(currLabelPair,'rows');
    Cset = unique(ic(:));
    countC = hist(ic(:),Cset);
    currCellOverlap = zeros(nanLabel,nanLabel);
    for j = 1:size(C,1)
        if j == 35
        end
        overlapEntry = C(j,:); 
        currCellOverlap(overlapEntry(1), overlapEntry(2)) = countC(j);
    end
    % in each of the matrices in the cellOverlap cell  array, the rows
    % indiate the previous cell ids and the columns indicate the current
    % cell ids
    cellOverlap{i} = currCellOverlap;
end


%% For each of the cellOverlap, get the cells which have changed ids
% each of the entries of cellChange is a cell array that keeps track of the
% changes along each of the label images
cellChange = cell(size(cellOverlap,1),1);
for i = 1:length(cellChange)
    currCellOverlap = cellOverlap{i};
    for row = 1:size(currCellOverlap)
        %See if any row has entries not equal to row numbers, meaning that
        %there are changed cell id's
        if( any(setdiff(find(currCellOverlap(row,:)), row) ) )
            ind = find(currCellOverlap(row,:));
            %cellChangeRecord is a struct that records the cell id in the
            %previous cellOverlap matrix and paired values that are the
            %cell ids in the current label + the pixels belonging to such cell id 
            cellChangeRecord.id = row;
            cellChangeRecord.record = zeros(length(ind),2);
            for z = 1:length(ind)
                couple = [ind(z) , currCellOverlap(row,ind(z))];
                cellChangeRecord.record(z,:) = couple;
            end
            if isempty(cellChange{i})
                cellChange{i} = {cellChangeRecord};
            else
                temp = cellChange{i};
                temp = {temp{:},cellChangeRecord};
                cellChange{i} = temp;
            end
        end
    end
end
%{
%% Find out the cells which change their pixel assignment more than 5 percent from the previous 
pixelChange = cell(size(cellChange,1),1);
for i =1:size(pixelChange,1)
    % loop through each of the structs in the current entry of celllChange 
    currCellChange = cellChange{i};
    for j = 1:length(currCellChange)
        currStruct = currCellChange{j};
        currID = currStruct.id;
        currRecord= currStruct.record;
        sumCurrRecord = sum(currRecord,1);
        % originalPixel is the original pixel size
        originalPixel = sumCurrRecord(2);
        % the row number of the cell id that the cell is supposed to stay
        % in 
        rowNum = find(currRecord(:,1)== currID);
        %changedPixels are the number of pixels that have switched
        changedPixels = originalPixel - currRecord(rowNum,2);
        % calculate the percentage of changed pixels. If more thatn 5
        % percent then add to pixelChange
        if(changedPixels/originalPixel > 0.05)
            temp = pixelChange{i};
            temp = [temp; [currID, changedPixels/originalPixel]];
            pixelChange{i} = temp;
        end
    end
end





%% take out the pairs that are the same
for i = 1:length(pairs)
  p = pairs{i};
  deletePair = zeros(1,1);
  for j = 1:size(p,1)
      if(p(j,1) ==p(j,2))
          deletePair(end,1) = j;
          deletePair(end+1,1) =0;
      end
  end
  deletePair = deletePair(1:end-1);
  p(deletePair,:) = [];
  pairs{i} = p;
end
%}
%{
%% Get the overlapping cells that are not nans
array = zeros(1,1);
for i = 1:length(pairs)
    p = pairs{i};
    p = [p(:,1);p(:,2)];
    array = [array;p];
    
    
end
array = unique(array);

%% Get the overlapping cells into groups
% cellList implements the graph data structure
cellList = cell(max(max(max(Lbl.Lbl))),1);
cellListEnd = 0;
% loop through pairs of overlap
for i = 1:length(pairs)
    currPairs = pairs{i};
    disp(strcat('i = ',num2str(i)));
    % for the current pair list, group cells into array
    for j = 1:size(currPairs,1)
        set = currPairs(j,:);
        % the rows to look up in currPairs
        rowsToLookup = zeros(size(currPairs,1),1);
        %rowsInd is the index of last nonzero element for rowsToLookup
        rowsInd = 0;
        % find out all the rows in currPairs that overlap with the variable
        % set
        for z = 1:size(currPairs,1)
            if any(intersect(set,currPairs(z,:)))
                rowsInd = rowsInd +1;
                rowsToLookup(rowsInd) = z;
            end   
        end
        %remove zero entries
        rowsToLookup = rowsToLookup(1:rowsInd);
        % form a new group of cells that overlap
        newGroup = unique([set,reshape(currPairs(rowsToLookup,:),1,[])]);
        % Check to see if the newGroup needs to be added to the existing
        % groups in cellList
        % the rows to look up in cellList
        rowsToLookup = zeros(size(cellList,1),1);
        %rowsInd is the index of last nonzero element for rowsToLookup
        rowsInd = 0;
        % find out all the rows in currPairs that overlap with the variable
        % set 
        for z = 1:cellListEnd
            if any(intersect(newGroup,cellList{z}))
                rowsInd = rowsInd +1;
                rowsToLookup(rowsInd) = z;
            end   
        end
        %remove zero entries
        rowsToLookup = rowsToLookup(1:rowsInd);
        % rowsToLookup ==0 is when the newGroup is unique and does not share overlap
        % with any existing group
        if ~any(find(rowsToLookup))
            cellListEnd = cellListEnd +1;
            cellList{cellListEnd} = newGroup;
        else
            % if the new group needs to be added to multiple groups, then combine those groups into one 
            if(length(rowsToLookup) >1)
                % sort the ib in ascending order
                rowsToLookup = sort(rowsToLookup);
                % combine the groups into one
                cellList{rowsToLookup(1)} = sort(unique([cellList{rowsToLookup(1)}, reshape(cell2mat(cellList(rowsToLookup(2:end))),1,[]),newGroup]));
                cellList{rowsToLookup(2:end)} = [];
                nonzeroCells = cellList(~cellfun('isempty',cellList));
                cellList= cell(size(cellList,1),1);
                cellList(1:size(nonzeroCells,1),1) = nonzeroCells;
                cellListEnd = size(nonzeroCells,1);
            else
                % add the new group of cell to a existing group
                cellList{rowsToLookup(1)} = unique([cellList{rowsToLookup(1)},newGroup]);
            end

        end
    end   
end
%}