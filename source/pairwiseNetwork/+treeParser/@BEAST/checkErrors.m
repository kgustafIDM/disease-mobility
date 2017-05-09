function checkErrors(obj)

for F = obj.fields,
    switch class(obj.(F{:}))
        case 'cell'
        otherwise
            error('BEAST:constructorRequiresCellInput','inputs must be of type cell')
    end
end

% special cases
    % leafNames
        if ~isempty(obj.leafNames) 
            if ~(length(obj.leafNames)==obj.numLeaves)
                error('BEAST:numleavesmismatch','number of leaf names must equal number of leaves')
            end
            if ~all(strcmp('char',cellfun(@class,obj.leafNames,'uniformoutput',false)))
                error('BEAST:leafNamesType','leaf names must be given as strings')
            end
        end

    % nodeNames
        if ~isempty(obj.nodeNames) 
            Len1=cellfun(@length,obj.nodeNames);
            if any(~(Len1==(2*obj.numLeaves-1))) 
                error('BEAST:numleavesmismatch','number of node names must equal number of nodes')
            end
        end

%         if ~isempty(obj.branchLengths)
%             if ~(length(obj.branchLengths)==obj.numTrees)
%                 error('BEAST:branchLengthInconsistent','number of branch length vectors does not match number of trees')
%             end
%             if ~all(cellfun(@length,obj.branchLengths)==(obj.numLeaves))
%                 error('BEAST:branchLengthsMissingData','number of branch lengths must match number of internal nodes + 1')
%             end
%         end
end