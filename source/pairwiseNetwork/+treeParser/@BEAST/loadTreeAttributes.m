function obj = loadTreeAttributes(obj,fidString,lastTipDate,lastTipDateUnits)
    % takes in BEAST .trees log file

    obj.fields={'fields','sourceFile','trees','isCoalescentTime','nodeTimes','numLineages',...
                'leafNames','leafNameTranslation','nodeNames','nodeAttributes'};
    
    obj.sourceFile=fidString;
    obj.timeUnits{1}=lastTipDateUnits;
    
%     [Missing, OutOfDate] = checkForMatlabToolbox('Bioinformatics Toolbox');

    % progress logicals
        step1=false;
        step2=false;
        step3=false;
        step4=false;
        step5=false;
        step6=false;

    fid = fopen(fidString,'r');
    numLines=fnumlines(fid);

    lineCounter=0;
    treeCounter=1;
    while lineCounter < numLines

        line=fgetl(fid);
        lineCounter=lineCounter+1;

        % finding number of leaves and pre-declaring node variables
            if ~isempty(regexp(line,'ntax=', 'once'))
                step1=true;
             
                NumLeaves = cell2mat(textscan(line,'%*s ntax=%d'));
                
                %pre-declare
                    obj.leafNames=cell(NumLeaves,1);
                    obj.leafNameTranslation=cell(NumLeaves,1); % leaves are labelled by integer in tree strings.  This is the lookup table to names

                ntaxaLine=lineCounter;
            end

        % leafNames
            if step1 && lineCounter > ntaxaLine + 1 && lineCounter <= (ntaxaLine + 1 + NumLeaves)
                tmp=textscan(line,'%s');
                obj.leafNames{lineCounter-ntaxaLine-1}=regexprep(tmp{1}{1},'''','');
                if lineCounter==(ntaxaLine + 1 + NumLeaves)
                    step2=true;
                end
            end

        % translation table
            if ~isempty(regexp(line,'Translate', 'once')),
                translateLine=lineCounter;
                step3=true;
            end
            if step1 && step2 && step3 && lineCounter > translateLine && lineCounter <= (translateLine + NumLeaves)
                tmp=textscan(line,'%n %[^,]');
                obj.leafNameTranslation{lineCounter-translateLine}={tmp{1},cell2mat(regexprep(tmp{2},'''',''))};
                if lineCounter==(translateLine+NumLeaves)
                    step4=true;
                end
            end

        % need bioinformatics toolbox to parse trees
%             if any(strcmp('Bioinformatics Toolbox',Missing))
%                 warning('bioinformatics toolbox missing. Cannot parse trees, so skipping.')
%             else
            
            % finding first tree in file
                if step1 && step2 && ~step5 && ~isempty(regexp(line,'tree ', 'once')) && isempty(regexp(line,'figtree', 'once')),
                    maxNumTrees=numLines-lineCounter;

                    % pre-declare
                        obj.nodeTimes=cell(maxNumTrees,1);
                        obj.nodeNames=cell(maxNumTrees,1);
                        obj.isCoalescentTime=cell(maxNumTrees,1);
                        obj.nodeAttributes=cell(maxNumTrees,1);

%                         if ~any(strcmp('Bioinformatics Toolbox',Missing))
%                             obj.trees=cell(maxNumTrees,1);
%                         end

                    step5=true;
                end

            % parsing trees
                if step1 && step2 && step5 && treeCounter <= maxNumTrees && ~isempty(regexp(line,'tree ', 'once')) && isempty(regexp(line,'figtree', 'once'))
                    
                    if treeCounter==1,
                        step6=true;
                    end
                    
                    % parsing tree string
                        nexusOut=cleanNexusTreeString(line);

                    % translating numbered nodes to names
                        if step4
                            for k=1:length(obj.leafNameTranslation)
                                nexusOut=regexprep(nexusOut,['(',num2str(obj.leafNameTranslation{k}{1}),'['],['(',num2str(obj.leafNameTranslation{k}{2}),'[']);
                                nexusOut=regexprep(nexusOut,[',',num2str(obj.leafNameTranslation{k}{1}),'['],[',',num2str(obj.leafNameTranslation{k}{2}),'[']);
                                nexusOut=regexprep(nexusOut,['(',num2str(obj.leafNameTranslation{k}{1}),':'],['(',num2str(obj.leafNameTranslation{k}{2}),':']);
                                nexusOut=regexprep(nexusOut,[',',num2str(obj.leafNameTranslation{k}{1}),':'],[',',num2str(obj.leafNameTranslation{k}{2}),':']);                            
                            end
                        end
                       
                    % import tree to matlab phytree object
                    v=ver;
                    if ~any(strcmpi('Bioinformatics Toolbox', {v.Name}))
                        error('requires bioinformatics toolbox to parse trees and build phytree object');
                    else
                        [tree, nodeAttributesOut] = phytreereadnexus(nexusOut);
                        obj.trees{treeCounter}=tree;    
                    end
                    % get number of leaves
                        numLeaves=get(tree,'NumLeaves');

                    % storing additional node attributes if they exist
                        if ~isempty(nodeAttributesOut)
                            obj.nodeAttributes{treeCounter}=nodeAttributesOut;
                        end                    
                        
                        % sometimes from saving imported info in FigTree, leaf attributes end up in the leaf names in []
                            if any(~cellfun(@isempty,regexp(obj.leafNames,'\[')))
                                names=get(tree,'nodenames');
                                fields={};
                                for c=1:length(names),
                                    nameIdx=find(~cellfun(@isempty,regexp(obj.leafNames,names{c})));
                                    if any(nameIdx)
                                        if numel(nameIdx)>1
                                            isCompleteName=false(length(nameIdx),1);
                                            for a=1:length(nameIdx)
                                                tmpSplit=get_tokens(obj.leafNames{nameIdx(a)},'\[');
                                                isCompleteName(a)=strcmp(names{c},tmpSplit{1});
                                            end
                                            nameIdx=nameIdx(isCompleteName);
                                            if numel(nameIdx)>1,
                                                error('leaf names are not unique.')
                                            end
                                        end
                                        nameBreakdown=get_tokens(obj.leafNames{nameIdx},'\[');
                                        tmpAtt=regexprep(regexprep(nameBreakdown{2}(1:end-1),',','&'),'"','');
                                        tmpBreakdown=get_tokens(tmpAtt,'&');
                                        for m=2:length(tmpBreakdown)
                                            tmp2=get_tokens(tmpBreakdown{m},'=');
                                            if c==1,
                                                if ~any(strcmp(tmp2{1},fieldnames(obj.nodeAttributes{treeCounter})))
                                                    fields{end+1}=tmp2{1};
                                                    obj.nodeAttributes{treeCounter}.(tmp2{1})=cell(2*numLeaves-1,1);
                                                else
                                                    warning(['node attribute, ',tmp2{1},' is defined in attributes and leaf name. ignoring leaf name specification']);
                                                end
                                            end
                                            obj.nodeAttributes{treeCounter}.(tmp2{1}){c}=tmp2{2};
                                        end
                                    end
                                end
                                for m=1:length(fields),
                                    emptyIdx=find(cellfun(@isempty,obj.nodeAttributes{treeCounter}.(fields{m})));
                                    if ~isempty(emptyIdx)
                                        for n=1:length(emptyIdx)
                                            obj.nodeAttributes{treeCounter}.(fields{m}){emptyIdx(n)}='';
                                        end
                                    end
                                end
                            end
                        
                        
                    % get node positions
                        phyData=phytreeNodePosition(tree);
                        tmpNodeTimes=[phyData.LeafDotsX; phyData.BranchDotsX];
                        tmpNodeNames=[get(tree,'LeafNames'); get(tree,'BranchNames')];
 
                        tmpIsCoalTime=[false(numLeaves,1);true(numLeaves-1,1)];

                    % align to known last tip date and cleanup
                        tmpNodeTimes=lastTipDate-max(tmpNodeTimes)+tmpNodeTimes;
                        tmpNodeTimes=datenumCleaner(tmpNodeTimes,lastTipDateUnits);

                    % sort dates
                        [tmpNodeTimes, nodeTimeSortIdx]=sort(tmpNodeTimes,'ascend');
                        obj.nodeTimes{treeCounter}=tmpNodeTimes;

                        tmpIsCoalTime = tmpIsCoalTime(nodeTimeSortIdx);
                        obj.isCoalescentTime{treeCounter} = tmpIsCoalTime;
                        obj.nodeNames{treeCounter} = tmpNodeNames(nodeTimeSortIdx);
                        
                        if ~isempty(nodeAttributesOut)
                            for F=reshape(fieldnames(obj.nodeAttributes{treeCounter}),1,[]),
                                obj.nodeAttributes{treeCounter}.(F{:})=obj.nodeAttributes{treeCounter}.(F{:})(nodeTimeSortIdx,:);
                            end
                        end
                        
                    treeCounter = treeCounter + 1;           
                    
                end

                % end of trees block
                    if step6 && ~isempty(regexp(line,'end;', 'once'))
                            obj.nodeTimes=obj.nodeTimes(1:(treeCounter-1));
                            obj.nodeNames=obj.nodeNames(1:(treeCounter-1));
                            obj.trees=obj.trees(1:(treeCounter-1));
                            obj.isCoalescentTime=obj.isCoalescentTime(1:(treeCounter-1));
                            obj.nodeAttributes=obj.nodeAttributes(1:(treeCounter-1));
                    end
                    
            end
%     end
    
    obj.getNumLineages;
    
    fclose(fid);
end
