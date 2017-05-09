function yearsOut = datenum2years(datenumIn)

switch class(datenumIn)
    case 'double'
        yearsOut=datenumIn/365.2431;
    case 'cell'
        yearsOut=datenumIn;
        for k=1:length(datenumIn)
            yearsOut{k}=datenumIn{k}/365.244;
        end
end
end