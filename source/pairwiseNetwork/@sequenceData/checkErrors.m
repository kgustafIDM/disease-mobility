function checkErrors(obj)
% sequenceData class
%

% type checking
for k=1:length(obj.cellFields),
    switch class(obj.(obj.cellFields{k}))
        case 'cell'
            % inside checking
%                 if ~isempty(obj.(obj.cellFields{k})) && ~all(strcmp('char',cellfun(@class,obj.(obj.cellFields{k}),'UniformOutput',false)))
%                     error('sequenceData:wrongDataType',['all entries in field ',obj.cellFields{k},' must be character strings'])
%                 end
                if k>1 && ~any([isempty(obj.(obj.cellFields{k-1})),isempty(obj.(obj.cellFields{k}))]) && ~all(size(obj.(obj.cellFields{k-1}))==size(obj.(obj.cellFields{k})))
                    error('sequenceData:mismatchedNumberOfEntries','number of entries must be the same for all fields')
                end
        otherwise
            error('sequenceData:wrongFieldType',['field type ',obj.cellFields{k},' must be cellarray'])
    end
end

for k=1:length(obj.doubleFields)
    switch class(obj.(obj.doubleFields{k}))
        case {'double','int32'}
            % check inside
                if ~isempty(obj.(obj.doubleFields{k})) && ~any(size(obj.(obj.doubleFields{k}))==1)
                    error('sequenceData:incorrectDimension',[obj.doubleFields{k},' must be single column vector'])
                end
                if ~isempty(obj.(obj.doubleFields{k})) && ~(length(obj.(obj.doubleFields{k}))==obj.numSeq)
                    error('sequenceData:mismatchedNumberOfEntries','number of entries must be the same for all fields')
                end
        otherwise
            error('sequenceData:wrongFieldType',['field type ',obj.doubleFields{k},' must be vector of doubles'])
    end
end
    

end