%Helper function to extract global decays from a masked region of a photon histogram.

function [decays] = maskRawData(rawData,masks)

%set up the variables needed
ADCres = size(rawData,3);
nDecays = size(masks,3);
decays = zeros(ADCres,nDecays);

%apply the masks to the images to get ROI-specific fluorescence decays
for j = 1:nDecays
    for i = 1:ADCres
        image = immultiply(masks(:,:,j),rawData(:,:,i));
        decays(i,j) = sum(sum(image));
    end
end

end

