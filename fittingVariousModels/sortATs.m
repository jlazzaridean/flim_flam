%This helper function ensures that taus are returned in ascending order.

function [aFs,tFs] = sortATs(aF,tF)

%transposing and shenanigans to set up an easily sortable array
sortable = zeros(size(aF,2),2);
sortable(:,1) = tF';
sortable(:,2) = aF';

%sort this combined table
sorted = sortrows(sortable);

%read out the values, still paired and now in order
aFs = sorted(:,2)';
tFs = sorted(:,1)';

end

