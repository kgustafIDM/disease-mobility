function datesOut = datenumCleaner(datesIn,dateUnits)
% takes approximate datenums in and converts them to integer datenums that
% can be compared with ==
%
% INPUT
%   datesIn: vector of dates
%   units: string either 'days' or 'years'
%

if nargin==1,
    dateUnits='days';
end
    
switch dateUnits
    case 'days'
        dateString=datestr(datesIn,'dd-mmm-yyyy');
    case 'years'
        dateString=datestr(datesIn*365.244,'dd-mmm-yyyy');
    otherwise
        error('datenumCleaner:IllegalUnits','units can only be days or years from 0 AD')
end

datesOut=datenum(dateString);


end