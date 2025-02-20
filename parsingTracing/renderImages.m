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
    % "photonsBin" will render a 32 bit photon image for analysis at the
    % binned spatial resolution
    % "photonsNoBin" will render the raw photon count from the tcscpc as a
    % 32 bit photon image for analysis at the native spatial resolution
    % "photonsScaleBin" will render an 8 bit photon image, scaled
    % for display, using the same binning used for fitting
    % "photonsScaleNoBin" will render an 8 bit photon image for display at
    % the native spatial resolution
    
%photonsMax = the maximum value to scale the photon count LUTs to in
%overlay images. If it's set to -1, it will use the max of the photons
%image.

function [] = renderImages(data,render,outPath,photonsMax,configS,varargin)

%varargin will be the range for each overlay as [limLow limHigh] in the
%order they are listed in the render matrix
nOverlays = nnz(contains(render,'overlay'));
if(nOverlays ~= size(varargin,2))
    error('Please enter a range as [low high] for each overlay requested in varargin');
end

%extract necessary params from the config file
ccvColor = configS(1,1).ccvColor;
model = configS(1,1).model;

for i=1:size(data,1)
    overlayCount = 1;
    for j=1:size(render,1)
        %extract the value to render & handle a few exception cases for
        %1exp (allow tau to be called a couple of different things)
        thisVal = char(render(j,1));
        if strcmp(model,'1exp')
            if strcmp(thisVal,'tau1') || strcmp(thisVal,'tm')
                thisVal = 'tau';
            end
        end
        
        %render the correct image based on the value entered for render
        oName = strcat(data(i,1).fName(1:end-4),'_',thisVal);
        if(contains(thisVal,'overlay'))
            ccvName = thisVal(9:end);
            if strcmp(model,'1exp') %handle the similar-sounding exception cases for 1exp fit.
                if strcmp(ccvName,'tau1') || strcmp(ccvName,'tm')
                    ccvName = 'tau';
                end
            end
            if(strcmp(data(i,1).binType,'std'))
                if(strcmp(ccvName,'tau1') || strcmp(ccvName,'tau2') || strcmp(ccvName,'tau3') || strcmp(ccvName,'tau4'))
                    tNum = str2double(ccvName(end));
                    [oFig, f2, oImg] = renderOverlay(data(i,1).photonsBinned,data(i,1).tau(:,:,tNum),photonsMax,varargin{1,overlayCount},ccvColor,oName);
                elseif (strcmp(ccvName,'a1') || strcmp(ccvName,'a2') || strcmp(ccvName,'a3') || strcmp(ccvName,'a4'))
                    aNum = str2double(ccvName(end));
                    [oFig, f2, oImg] = renderOverlay(data(i,1).photonsBinned,data(i,1).a(:,:,aNum),photonsMax,varargin{1,overlayCount},ccvColor,oName);
                else
                    [oFig, f2, oImg] = renderOverlay(data(i,1).photonsBinned,data(i,1).(ccvName),photonsMax,varargin{1,overlayCount},ccvColor,oName);
                end
            elseif(strcmp(data(i,1).binType,'bh'))
                if(strcmp(ccvName,'tau1') || strcmp(ccvName,'tau2') || strcmp(ccvName,'tau3') || strcmp(ccvName,'tau4'))
                    tNum = str2double(ccvName(end));
                    [oFig, f2, oImg] = renderOverlay(data(i,1).photons,data(i,1).tau(:,:,tNum),photonsMax,varargin{1,overlayCount},ccvColor,oName);
                elseif (strcmp(ccvName,'a1') || strcmp(ccvName,'a2') || strcmp(ccvName,'a3') || strcmp(ccvName,'a4'))
                    aNum = str2double(ccvName(end));
                    [oFig, f2, oImg] = renderOverlay(data(i,1).photons,data(i,1).a(:,:,aNum),photonsMax,varargin{1,overlayCount},ccvColor,oName);
                else
                    [oFig, f2, oImg] = renderOverlay(data(i,1).photons,data(i,1).(ccvName),photonsMax,varargin{1,overlayCount},ccvColor,oName);
                end                
            else
                error('Unrecognized bin type. Cannot render overlay');
            end
            saveas(oFig,strcat(outPath,filesep,oName,'.pdf'),'pdf');
            imwrite(oImg.cdata,strcat(outPath,filesep,oName,'.tiff'));
            close(f2);
            overlayCount = overlayCount + 1;
        elseif(strcmp(thisVal,'photonsBin'))
            writeTIFF(data(i,1).photonsBinned,'analysis',outPath,oName);
        elseif(strcmp(thisVal,'photonsNoBin'))
            writeTIFF(data(i,1).photons,'analysis',outPath,oName);
        elseif(strcmp(thisVal,'photonsScaleNoBin'))
            %scale the unbinned photons image and write to 8 bit TIFF
            if(photonsMax == -1)
                photonsMax = max(data(i,1).photons,[],'all');
            end
            scaled = data(i,1).photons ./ photonsMax .* 255;
            writeTIFF(scaled,'photonsScaled',outPath,oName);
        elseif(strcmp(thisVal,'photonsScaleBin'))
            %scale the binned photons image and write to 8 bit TIFF
            if(photonsMax == -1)
                photonsMax = max(data(i,1).photonsBinned,[],'all');
            end
            scaled = data(i,1).photonsBinned ./ photonsMax .* 255;
            writeTIFF(scaled,'photonsScaled',outPath,oName);           
        elseif(strcmp(thisVal,'tau1') || strcmp(thisVal,'tau2') || strcmp(thisVal,'tau3') || strcmp(thisVal,'tau4'))
            tNum = str2double(thisVal(end));
            writeTIFF(data(i,1).tau(:,:,tNum),'analysis',outPath,oName);
        elseif(strcmp(thisVal,'a1') || strcmp(thisVal,'a2') || strcmp(thisVal,'a3') || strcmp(thisVal,'a4'))
            aNum = str2double(thisVal(end));
            writeTIFF(data(i,1).a(:,:,aNum),'analysis',outPath,oName);
        else
            %if the field hasn't been handled yet, it msut be an exact
            %match or else the code will throw and error.
            if isfield(data,thisVal)
                writeTIFF(data(i,1).(thisVal),'analysis',outPath,oName);
            else
                error(strcat("Data structure does not contain ", thisVal, " - cannot be rendered."));
            end
        end
    end
end

end

