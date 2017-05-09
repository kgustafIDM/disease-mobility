function [nexusStringOut] = cleanNexusTreeString(beastTreeLineIn)
% strips nexus tree line string from BEAST .trees file down to a parsable nexus string
% DEPENDENCIES: get_tokens and get_cell_tokens from fileExchange
%


firstStringBreakdown = get_tokens(beastTreeLineIn,'\[');

secondStringBreakdown = get_cell_tokens(firstStringBreakdown,'\]');

% finding opening parenthesis
    firstTreeBin=0; ii=0;
    while firstTreeBin==0,
        ii=ii+1;
        tmp = regexp(secondStringBreakdown{ii},'(', 'once');
        if any(~cellfun(@isempty,tmp))
            cellBin=find(~cellfun(@isempty,tmp),1);
            firstTreeBin=ii;
        end
    end

% build newick string
    nexusStringOut = [];
    for k=firstTreeBin:(length(secondStringBreakdown)-1),
        nexusStringOut = [nexusStringOut, secondStringBreakdown{k}{cellBin},'[',regexprep(secondStringBreakdown{k+1}{cellBin-1},',','&'),']'];
    end
    nexusStringOut = [nexusStringOut, secondStringBreakdown{end}{2}];
    
end