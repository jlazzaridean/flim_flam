%Writes the array specified by the input parameter data to TIFF. 

%The file is saved to filename within path.

function [] = writeTIFF(data,type,path,filename)

if(size(data,3) >1)
    error('You are trying to render a parameter with 3 dimensions, likely a or tau rather than an individual component.');
end

%reflect the data on the y axis because matlab flips images for some reason
data1 = flip(data,1);

fs=filesep;
outputFileName = strcat(path,fs,filename,'.tiff');
t = Tiff(outputFileName,'w');
% Setup tags
% Lots of info here:
% https://www.mathworks.com/help/matlab/ref/tiff.html
tagstruct.ImageLength     = size(data1,1);
tagstruct.ImageWidth      = size(data1,2);
tagstruct.RowsPerStrip    = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software        = 'MATLAB';

%set up the data formatting either for analysis or display
if(strcmp(type,'photonsScaled')) %save an 8 bit TIFF on a fixed scale
    tagstruct.SampleFormat = 1; %uint data
    tagstruct.Photometric = 1; %grayscale, min is black
    tagstruct.BitsPerSample   = 8; %8 bit TIFF
    tagstruct.SamplesPerPixel = 1; %there is one dimension in the scaled photons
    data1 = uint8(data1);
elseif(strcmp(type,'analysis')) %retain all of the information for analysis
    tagstruct.SampleFormat = 3; %IEEE floating point
    tagstruct.Photometric     = 1; %grayscale min is black
    tagstruct.BitsPerSample   = 32; %32 bit TIFF
    tagstruct.SamplesPerPixel = 1; %there is 1 dimension
    data1 = single(data1);
end


t.setTag(tagstruct);
t.write(data1);
t.close();

end

