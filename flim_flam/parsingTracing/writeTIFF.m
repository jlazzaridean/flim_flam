%Writes the array specified by the input parameter data to TIFF. 

%The file is saved to filename within path.

function [] = writeTIFF(data,path,filename)
fs=filesep;
outputFileName = strcat(path,fs,filename,'.tiff');
t = Tiff(outputFileName,'w');
% Setup tags
% Lots of info here:
% https://www.mathworks.com/help/matlab/ref/tiff.html
tagstruct.ImageLength     = size(data,1);
tagstruct.ImageWidth      = size(data,2);
tagstruct.Photometric     = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample   = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';
tagstruct.SampleFormat = 3;
t.setTag(tagstruct);
t.write(single(data));
t.close();

end

