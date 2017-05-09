function checkErrors(obj)

% type checking
    switch class(obj.timeUnits)
        case 'char'
            if ~(strcmp('years',obj.timeUnits) || strcmp('days',obj.timeUnits)),
                error('baseClass:timeUnitsType','time units must be ''years'' or ''days'' only')
            end
        otherwise
            error('baseClass:timeUnits','time units must be character string')
    end
    switch class(obj.nodeTimes)
        case 'cell'
        otherwise
            error('baseClass:constructorRequiresCellInput','inputs must be a cellarray of vectors')
    end
    switch class(obj.isCoalescentTime)
        case 'cell'
            for k=1:length(obj.isCoalescentTime),
                switch class(obj.isCoalescentTime{k})
                    case 'logical'
                    otherwise
                        error('baseClass:isCoalescentTimeMustBeLogical','isCoalescentTime must contain logical vectors')
                end
            end
        otherwise
            error('baseClass:constructorRequiresCellInput','inputs must be a cellarray of vectors')                
    end

% data match checking
    if any(~(size(obj.nodeTimes)==size(obj.isCoalescentTime))),
        error('baseClass:lengthsMustMatch','nodeTimes and isCoalescentTime must have same number of cells')
    end
    % inside
        Len1=cellfun(@length,obj.nodeTimes);
        Len2=cellfun(@length,obj.isCoalescentTime);
        if any(~(Len1==Len2)),
            error('basicClass:constructorLengthMismatch','all vectors at a cell position must have matching lengths')
        end

end