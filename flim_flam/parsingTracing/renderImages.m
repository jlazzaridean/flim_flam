%This function renders various arrays from pixel to pixel fit data as images and saves them as TIFF files to outPath.

%The data input to this function is the result of pxwiseAnalysis. This is
%the same function called within pxwiseAnalysis.
%    The function looks for the string binType to correctly process the
%    overlay

%The parameters to render are contained within the N x 1 array of strings.
%They must be either direct matches to fields in the pixelwiseAnalysis or
%in the format "overlay_field" where field is a field in pixelwiseAnalysis.
%   For individual components of lifetime or coefficients, specify them as
%   tau# or a#, e.g. tau1 or a2.

%For each overlay that is requested, an additional varargin (separated by
%commas) must be used to specify the desired range. Each input argument
%should be in the format [limLow limHigh].

%Bin type is 'bh' or 'std'. This must be specified for overlays to render
%correctly. 
%   'bh' Becker-Hickl style upsampling/moving average binning, lifetime
%   image will be display over the unbinned photon image
%   'std' conventional binning where bin 2 makes 2x2 sets of px into 1 px.
%   lifetime immage will be displayed over the binned photon image.

%The render input is an array of strings (dimensions Nx1) with the names of
%the variables to render. These are fields in the data structure. For
%different tau components, use "tau1", "tau2", etc. For different
%coefficients, use "a1", "a2", etc. All other values should be called
%directly as their names.
    % "photons" will render the photon count on the binned image
    % "photons_raw" will render the raw photon count from the tcscpc
    
%photonsMax = the maximum value to scale the photon count LUTs to in
%overlay images. If it's set to -1, it will use the max of the photons
%image.

function [] = renderImages(data,render,outPath,photonsMax,varargin)

%varargin will be the range for each overlay as [limLow limHigh] in the
%order they are listed in the render matrix
nOverlays = nnz(contains(render,'overlay'));
if(nOverlays ~= size(varargin,2))
    error('Please enter a range as [low high] for each overlay requested in varargin');
end

for i=1:size(data,1)
    overlayCount = 1;
    for j=1:size(render,1)
        thisVal = char(render(j,1));
        oName = strcat(data(i,1).fName(1:end-4),'_',thisVal);
        if(contains(thisVal,'overlay'))
            ccvName = thisVal(9:end);
            if(strcmp(data(i,1).binType,'std'))
                if(strcmp(ccvName,'tau1') || strcmp(ccvName,'tau2') || strcmp(ccvName,'tau3') || strcmp(ccvName,'tau4'))
                    tNum = str2double(ccvName(end));
                    oFig = renderOverlay(data(i,1).photonsBinned,data(i,1).tau(:,:,tNum),photonsMax,varargin{1,overlayCount},oName);
                elseif (strcmp(ccvName,'a1') || strcmp(ccvName,'a2') || strcmp(ccvName,'a3') || strcmp(ccvName,'a4'))
                    aNum = str2double(ccvName(end));
                    oFig = renderOverlay(data(i,1).photonsBinned,data(i,1).a(:,:,aNum),photonsMax,varargin{1,overlayCount},oName);
                else
                    oFig = renderOverlay(data(i,1).photonsBinned,data(i,1).(ccvName),photonsMax,varargin{1,overlayCount},oName);
                end
            elseif(strcmp(data(i,1).binType,'bh'))
                if(strcmp(ccvName,'tau1') || strcmp(ccvName,'tau2') || strcmp(ccvName,'tau3') || strcmp(ccvName,'tau4'))
                    tNum = str2double(ccvName(end));
                    oFig = renderOverlay(data(i,1).photons,data(i,1).tau(:,:,tNum),photonsMax,varargin{1,overlayCount},oName);
                elseif (strcmp(ccvName,'a1') || strcmp(ccvName,'a2') || strcmp(ccvName,'a3') || strcmp(ccvName,'a4'))
                    aNum = str2double(ccvName(end));
                    oFig = renderOverlay(data(i,1).photons,data(i,1).a(:,:,aNum),photonsMax,varargin{1,overlayCount},oName);
                else
                    oFig = renderOverlay(data(i,1).photons,data(i,1).(ccvName),photonsMax,varargin{1,overlayCount},oName);
                end                
            else
                error('Unrecognized bin type. Cannot render overlay');
            end
            saveas(oFig,strcat(outPath,filesep,oName,'.pdf'),'pdf');
            overlayCount = overlayCount + 1;
        elseif(strcmp(thisVal,'photons'))
            writeTIFF(data(i,1).photonsBinned,outPath,oName);
        elseif(strcmp(thisVal,'photons_raw'))
            writeTIFF(data(i,1).photons,outPath,oName);            
        elseif(strcmp(thisVal,'tau1') || strcmp(thisVal,'tau2') || strcmp(thisVal,'tau3') || strcmp(thisVal,'tau4'))
            tNum = str2double(thisVal(end));
            writeTIFF(data(i,1).tau(:,:,tNum),outPath,oName);
        elseif(strcmp(thisVal,'a1') || strcmp(thisVal,'a2') || strcmp(thisVal,'a3') || strcmp(thisVal,'a4'))
            aNum = str2double(thisVal(end));
            writeTIFF(data(i,1).a(:,:,aNum),outPath,oName);
        else
            writeTIFF(data(i,1).(thisVal),outPath,oName);
        end
    end
end

end

