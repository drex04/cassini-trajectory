function [JD,MJD] = DateConvert(year,month,day,hour,minute,second)
% Author: Drew Nollsch
% Date: 9/8/14
% ASTE 580, Professor Michael Gabor
%
% A function to calculate Julian Date and Modified Julian Date from
% calendar date and time
%
% Inputs:   year - calendar year
%           month - calendar month (numeral or string)
%           day - calendar day
%           hour
%           minute
%           second
%
% Outputs:  JD - Julian Date, days
%           MJD - Modified Julian Date, days
%%
% Allow for string inputs for month
if strcmpi(month,'January')
    month = 1;
elseif strcmpi(month,'February')
    month = 2;
elseif strcmpi(month,'March')
    month = 3;
elseif strcmpi(month,'April')
    month = 4;
elseif strcmpi(month,'May')
    month = 5;
elseif strcmpi(month,'June')
    month = 6;
elseif strcmpi(month,'July')
    month = 7;
elseif strcmpi(month,'August')
    month = 8;
elseif strcmpi(month,'September')
    month = 9;
elseif strcmpi(month,'October')
    month = 10;
elseif strcmpi(month,'November')
    month = 11;
elseif strcmpi(month,'December')
    month = 12; 
end

% Compute Julian Date
JD = 367*year - floor(7*(year+floor((month+9)/12))/4) + floor(275*month/9) + day + 1721013.5 + ((second/60+minute)/60+hour)/24;

% Compute Modified Julian Date
MJD = JD - 2400000.5;
