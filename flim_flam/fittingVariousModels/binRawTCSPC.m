%Bins an inputImg in both space and time, with options for different bin styles.

%The function returns a binnedImg and a binnedIRF (binned version of the
%input IRF, which must be specified for consistent processing of the time
%dimension during exponential fitting).

%inputImg is an X by Y by N array. inputIRF is an N by 1 array. 

%spatialBin is the extent of spatial binning; what this means depends on
%binType. binType is either 'bh' or 'std'.
%   'bh' - Becker and Hickl style binning as implemented in SPCImage;
%   resulting pixels are the sum of all pixels that are spatialBin away.
%   Therefore, bin 1 means that each pixel is a moving sum of its value
%   plus the eight surrounding pixels, but the binned image is the same
%   size as the original.
%   'std' - Conventional microscopy binning (also binning as implemented in
%   SymPhoTime). Resulting pixels are the sum spatialBin x spatialBin
%   pixels, and the resulting image is smaller by a factor of spatialBin.

%The N dimension (ADC) will be binned by adcBin in a standard way
%(summation of adjacent time bins). If ADC binning is used the IRF must be
%binned to match the image.

%Written by Julia Lazzari-Dean, Last edits 10/12/2019.

function [binnedImg,binnedIRF] = binRawTCSPC(inputImg,inputIRF,spatialBin,binType,adcBin)

startR = size(inputImg,1);
startC = size(inputImg,2);
startADC = size(inputImg,3);

spatialBin = round(spatialBin); % round the binning if non integer
adcBin = round(adcBin);
finalADC = floor(startADC/adcBin);

if (strcmp(binType,'std'))
    finalR = floor(startR/spatialBin); %if finalR is non integer, round down
    finalC = floor(startC/spatialBin);
    binnedImg = zeros(finalR,finalC,finalADC);
    for i = 1:finalR
        for j = 1:finalC
            if adcBin == 1
                values = inputImg(i*spatialBin-spatialBin+1:i*spatialBin,j*spatialBin-spatialBin+1:j*spatialBin,:);
                binnedImg(i,j,:) = sum(values,[1 2]);
            elseif adcBin > 1
                for k = 1:finalADC
                    values = inputImg(i*spatialBin-spatialBin+1:i*spatialBin,j*spatialBin-spatialBin+1:j*spatialBin,...
                        adcBin*k-adcBin+1:adcBin*k);
                    binnedImg(i,j,k) = sum(values,'all');
                end
            else
                error('ADC bin must be an integer value greater than or equal to 1.');
            end
        end
    end
elseif (strcmp(binType,'bh'))
    binnedImg = zeros(startR,startC,finalADC);
    for i = 1:startR
        for j = 1:startC
            %determine the boundaries of the area to sum
            minR = i - spatialBin;
            if(minR < 1)
                minR = 1;
            end
            minC = j - spatialBin;
            if(minC < 1)
                minC = 1;
            end
            maxR = i + spatialBin;
            if(maxR > startR)
                maxR = startR;
            end
            maxC = j + spatialBin;
            if(maxC > startC)
                maxC = startC;
            end
            if adcBin == 1
                values = inputImg(minR:maxR,minC:maxC,:);
                binnedImg(i,j,:) = sum(values,[1 2]);
            elseif adcBin > 1
                for k = 1:finalADC
                    values = inputImg(minR:maxR,minC:maxC,adcBin*k-adcBin+1:adcBin*k);
                    binnedImg(i,j,k) = sum(values,'all');
                end
            else
                error('ADC bin must be an integer value greater than or equal to 1.');
            end
        end
    end
end

%ADC binning only makes sense as 'standard' (no spatial moving average)
if (adcBin > 1)
    %bin the IRF in ADC space
    binnedIRF = zeros(finalADC,1);
    for i = 1:finalADC
        binnedIRF(i,1) = sum(inputIRF(adcBin*i-adcBin+1:adcBin*i,1));
    end
else
    binnedIRF = inputIRF;
end

end

