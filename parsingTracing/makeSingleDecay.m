%Helper function to extract a single global decay from all pixels of a photon histogram.

function [decay] = makeSingleDecay(tcspc,offset)

decay = -1;
for j = 1:size(tcspc,3)
    decay(j,1) = sum(sum(tcspc(:,:,j)));
end

%subtract the offset (or zero) from each position
decay = decay - offset;
decay (decay < 0) = 0;

end

