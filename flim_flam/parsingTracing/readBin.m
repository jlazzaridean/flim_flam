% Reads PicoQuant bin files exported from SymPhoTime64
% This is demo code. Use at your own risk.
% No warranties and no support.
% Peter Kapusta, PicoQuant GmbH, January 22th, 2018

%adapted/truncated by JLD 10/4/19

function [dataHist] = readBin(fullFileName)

fid=fopen(fullFileName);

PixX = fread(fid, 1, 'int32'); %x = number of COLUMNS
PixY = fread(fid, 1, 'int32'); %y = number of ROWS
PixResol = fread(fid, 1, 'float32');
TCSPCChannels = fread(fid, 1, 'int32');
TimeResol = fread(fid, 1, 'float32');

dataHist=zeros(PixY,PixX,TCSPCChannels);

for y=1:PixY
    for x=1:PixX
        dataHist(y,x,:)=fread(fid, TCSPCChannels, 'int32');
    end
end

fclose(fid);

end

