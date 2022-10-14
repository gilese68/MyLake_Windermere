% PhysPar for lake_params
% Import from Excel - need to move this file somewhere else
C = readcell('lake_params.xlsx');

lake_params = cell(61,2);
lake_params(:,1) = {cell(1,(size(C,2)-1))};

for i = 1:(size(C,2)-1)
    for ii = 1:size(C,1)
        lake_params{ii, 1}{1, i} = C{ii, i+1};
    end
end

% Import parameter names
for i = 1:size(C,1)
lake_params{i,2} = C{i,1};
end