%Reads a .csv config file containing fit settings and acquisition parameters

%The optional input argument can be used to specify the directory opened by
%uigetfile.

%The function returns a 1x1 structure containing the information in the
%configuration file and the name of the configuration file so that the user
%can approve it later.

function [config, fName] = readConfig(varargin)

%open the convenient path if it is provided to the function
fsep = filesep;
if(size(varargin,2) == 1)
    pathToOpen = varargin{1,1};
else
    pathToOpen = pwd;
end
if(pathToOpen(end) ~= fsep)
    pathToOpen(end+1) = fsep;
end

%prompt the user for a config file .csv (basically TCSPC system parameters).
%read this information into a structure that can be parsed by the other
%scripts so that all variables don't have to be entered individually
[fName,pName] = uigetfile([pathToOpen '*.csv'],'Please select .csv CONFIG file');
fullPath = [pName fsep fName];

T = readtable(fullPath);
configLong = table2struct(T);

%make sure that the correct fields are present in the config file
fields = fieldnames(configLong);
paramFields = contains(fields,'Param');

switch configLong(1,1).model
    case '1exp'
        fixed = -1;
        start = configLong(1,1).startParam;
    case '2exp'
        fixed = configLong(1,1).fixedParam_1;
        fixed(1,2) = configLong(1,1).fixedParam_2;
        fixed(1,3) = configLong(1,1).fixedParam_3;
        fixed(1,4) = configLong(1,1).fixedParam_4;
        start = configLong(1,1).startParam_1;
        start(1,2) = configLong(1,1).startParam_2;
        start(1,3) = configLong(1,1).startParam_3;
        start(1,4) = configLong(1,1).startParam_4;
    case '3exp'
        fixed = configLong(1,1).fixedParam_1;
        fixed(1,2) = configLong(1,1).fixedParam_2;
        fixed(1,3) = configLong(1,1).fixedParam_3;
        fixed(1,4) = configLong(1,1).fixedParam_4;
        fixed(1,5) = configLong(1,1).fixedParam_5;
        fixed(1,6) = configLong(1,1).fixedParam_6;
        start = configLong(1,1).startParam_1;
        start(1,2) = configLong(1,1).startParam_2;
        start(1,3) = configLong(1,1).startParam_3;
        start(1,4) = configLong(1,1).startParam_4;
        start(1,5) = configLong(1,1).startParam_5;
        start(1,6) = configLong(1,1).startParam_6;
    case '4exp'
        fixed = configLong(1,1).fixedParam_1;
        fixed(1,2) = configLong(1,1).fixedParam_2;
        fixed(1,3) = configLong(1,1).fixedParam_3;
        fixed(1,4) = configLong(1,1).fixedParam_4;
        fixed(1,5) = configLong(1,1).fixedParam_5;
        fixed(1,6) = configLong(1,1).fixedParam_6;
        fixed(1,7) = configLong(1,1).fixedParam_7;
        fixed(1,8) = configLong(1,1).fixedParam_8;
        start = configLong(1,1).startParam_1;
        start(1,2) = configLong(1,1).startParam_2;
        start(1,3) = configLong(1,1).startParam_3;
        start(1,4) = configLong(1,1).startParam_4;
        start(1,5) = configLong(1,1).startParam_5;
        start(1,6) = configLong(1,1).startParam_6;
        start(1,7) = configLong(1,1).startParam_7;
        start(1,8) = configLong(1,1).startParam_8;
    otherwise
        error('Unrecognized model field');
end


config = rmfield(configLong,fields(paramFields));

config(1,1).stFi = [config(1,1).fit_start config(1,1).fit_end];
config(1,1).IRFlim = [config(1,1).irf_start config(1,1).irf_end];
config = rmfield(config,{'fit_start','fit_end','irf_start','irf_end'});

config(1,1).startParam = start;
config(1,1).fixedParam = fixed;

end

