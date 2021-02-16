%Bins an input global decay and IRF in the ADC dimension.

%The function returns a list of binnedDecays and a binnedIRF (binned version of the
%input IRF, which must be specified for consistent processing of the time
%dimension during exponential fitting).

%inputDecays is an N by X array of X decays with N ADC time bins where each row is an ADC time bin and 
%each column is a unique decay. inputIRF is an N by 1 array. 

%The N dimension (ADC) will be binned by adcBin in a standard way
%(summation of adjacent time bins). If ADC binning is used the IRF must be
%binned to match the image.

%Written by Julia Lazzari-Dean, Last edits 10/14/2019.

function [binnedDecays,binnedIRF] = binADCGlobal(inputDecays,inputIRF,adcBin)

startADC = size(inputDecays,1);
adcBin = round(adcBin);
finalADC = floor(startADC/adcBin);
nDecays = size(inputDecays,2);
binnedDecays = zeros(finalADC,nDecays);
binnedIRF = zeros(finalADC,1);

for i = 1:nDecays
    %ADC binning only makes sense as 'standard' (no spatial moving average)
    if (adcBin > 1)
        %bin the IRF in ADC space
        for j = 1:finalADC
            binnedIRF(j,1) = sum(inputIRF(adcBin*j-adcBin+1:adcBin*j,1));
            binnedDecays(j,i) = sum(inputDecays(adcBin*j-adcBin+1:adcBin*j,i));
        end
    else
        binnedIRF = inputIRF;
        binnedDecays = inputDecays;
    end
end

end

