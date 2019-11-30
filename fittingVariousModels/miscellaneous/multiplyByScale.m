function [scaledDecays] = multiplyByScale(avgDecays,roiScales)

%check dimensions - roiScales should be a single column (as exported from ImageJ)
%decays should be a series of columns (1 column matches a row of the sizes)
nDecays = size(avgDecays,2);
nSizes = size(roiScales,1);

if(nDecays ~= nSizes)
    error('Dimension mismatch.');
end

scaledDecays = zeros(size(avgDecays));
ADCres = size(avgDecays,1);

for i = 1:nDecays
    for j = 1:ADCres
        scaledDecays(j,i) = avgDecays(j,i)*roiScales(i,1);
    end    
end


end

