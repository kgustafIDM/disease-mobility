function exportBEASTXML(obj,fileId,varargin)
% inputs: fileId with extension (.xml unless you've got a reason to do otherwise)
%
%   varargin :: key-value pairs
%       'CaseID' :: 'FastaName' (default) or 'EPID'
%       'isolateDateUnits' :: 'days' (default) or 'years'
%       'locationIdType' :: 'name' (default)  or 'shapeFileId'
%

% vararing handling
    options=struct('CaseID','FastaName','isolateDateUnits','days','locationIdType','name','copyModelConfig',false,'overwrite',false);
    inputs=varargin;
    options = keyValuePairVararginHandler(options,inputs);


s=struct();

% source data files
    s.beast.Comment{1}.Text = ['Source Files:: ',obj.sourceFile];

% including constructorQuery
    comment='Constructor Database Query:: ';
    for k=1:length(obj.constructorQuery)
            switch class(obj.constructorQuery{k})
                case 'cell'
                    for m=1:length(obj.constructorQuery{k})
                        switch class(obj.constructorQuery{k}{m})
                            %case 'cell'
                            case 'char'
                                tmp = ['_',obj.constructorQuery{k}{m}];
                            case 'double'
                                tmp = regexprep(num2str(obj.constructorQuery{k}{m},'_%d'),' ','');
                        end
                    end
                case 'char'
                    tmp = ['_',obj.constructorQuery{k}];
                case 'double'
                    tmp = regexprep(num2str(obj.constructorQuery{k},'_%d'),' ','');
            end
            comment = [comment tmp];
    end
    s.beast.Comment{2}.Text = comment;
    s.beast.Comment{3}.Text = ['This file was created on:: ',datestr(now,'dd-mmm-yyyy HH:MM:SS')];

s.beast.taxa.Attributes.id='taxa';
s.beast.alignment.Attributes.id='alignment';
s.beast.alignment.Attributes.dataType='nucleotide';

if ~iscell(obj.(options.CaseID)) || ~ischar(obj.(options.CaseID){1}) 
    tmp=cell(length(obj.(options.CaseID)),1);
    for k=1:length(obj.(options.CaseID))
        tmp{k}=obj.(options.CaseID)(k,:);
    end
    obj.(options.CaseID)=tmp;
end
    
for k=1:obj.numSeq,
    
    % taxa block
        s.beast.taxa.taxon{k}.Attributes.id = regexprep(obj.(options.CaseID){k},'-','_');
        switch options.isolateDateUnits
            case 'days'
                s.beast.taxa.taxon{k}.date.Attributes.value=num2str(obj.isolateDate(k));
            case 'years'
                s.beast.taxa.taxon{k}.date.Attributes.value=num2str(datenum2years(obj.isolateDate(k)));
        end
        s.beast.taxa.taxon{k}.date.Attributes.direction='forwards';
        s.beast.taxa.taxon{k}.date.Attributes.units=options.isolateDateUnits;
        
        % optional data
            % location
                AdminLevelCounter=1;
                if ~isempty(obj.AdminL0)
                    s.beast.taxa.taxon{k}.attr{AdminLevelCounter}.Attributes.name='AdminL0';
                    s.beast.taxa.taxon{k}.attr{AdminLevelCounter}.Text = lower(obj.AdminL0{k});
                    AdminLevelCounter=AdminLevelCounter+1;
                end
                   
                if ~isempty(obj.AdminL1)
                    s.beast.taxa.taxon{k}.attr{AdminLevelCounter}.Attributes.name='AdminL1';
                    s.beast.taxa.taxon{k}.attr{AdminLevelCounter}.Text = lower(obj.AdminL1{k});
                    AdminLevelCounter=AdminLevelCounter+1;
                end
                   
                if ~isempty(obj.AdminL2)
                    s.beast.taxa.taxon{k}.attr{AdminLevelCounter}.Attributes.name='AdminL2';
                    s.beast.taxa.taxon{k}.attr{AdminLevelCounter}.Text = lower(obj.AdminL2{k});
                    AdminLevelCounter=AdminLevelCounter+1;
                end
                
                if ~isempty(obj.ClusterID)
                    s.beast.taxa.taxon{k}.attr{AdminLevelCounter}.Attributes.name='ClusterID';
                    s.beast.taxa.taxon{k}.attr{AdminLevelCounter}.Text = lower(obj.ClusterID{k});
                end
    
    % alignment block
        s.beast.alignment.sequence{k}.taxon.Attributes.idref = s.beast.taxa.taxon{k}.Attributes.id;
        s.beast.alignment.sequence{k}.Text = upper(obj.sequence{k});
    
    
    
end

% parse filename
    [fileDir,fileName] = fileparts(fileId);
            
% copy model config
    if ~islogical(options.copyModelConfig) && ischar(options.copyModelConfig)
        % provenence
            s.beast.Comment{1}.Text = [s.beast.Comment{1}.Text,'; ',options.copyModelConfig];

        % import model config from BEAST xml
            modelConfig =  xml2struct(options.copyModelConfig);
            if ~isfield(modelConfig,'beast')
                error('model config file must be BEAST xml')
            end

        % append to output struct
            fields=fieldnames(modelConfig.beast);
            for k=1:length(fields),
                if ~strcmp('taxa',fields{k}) && ~strcmp('alignment',fields{k})
                    s.beast.(fields{k})=modelConfig.beast.(fields{k});
                end                
            end
        
        % replace file names in imported config with current file name
        
            % are log files with the same name already in the directory?
            % if yes, iterate filename.
            if ~options.overwrite
                if isdir(fileDir)
                    tmp=dir(fileDir);
                    tmpIdx=find(cellfun(@any,regexp({tmp.name},'.log')) | cellfun(@any,regexp({tmp.name},'.trees')) | cellfun(@any,regexp({tmp.name},'.ops')));
                    if ~isempty(tmpIdx)
                        tmpName={};
                        for k=1:length(tmpIdx)
                            [~, fN] = fileparts(tmp(tmpIdx(k)).name);
                            tmpName{end+1}=fN;
                        end
                        nameRecursion=true;
                        newFileName=fileName;
                        newNumber=1;
                        while nameRecursion
                            if any(strcmp(newFileName,tmpName))
                                newFileName=[fileName,'_',num2str(newNumber)];
                                newNumber=newNumber+1;
                            else
                                nameRecursion=false;
                                fileName=newFileName;
                                fileId=[fileDir,'\',fileName];
                            end
                        end
                    end
                end
            end
                
            % insert fileName
                s.beast.mcmc.log{2}.Attributes.fileName=[fileName,'.log'];
                s.beast.mcmc.logTree.Attributes.fileName=[fileName,'.trees'];
            
                % are we creating an operator Analysis file?
                    if isfield(s.beast.mcmc.Attributes,'operatorAnalysis')
                        s.beast.mcmc.Attributes.operatorAnalysis=[fileName,'.ops'];
                    end
    elseif islogical(options.copyModelConfig) && ~options.copyModelConfig
    else
        error('copyModelConfig option must be a filename string')
    end



% write file
    if ~isdir(fileDir) && ~isempty(fileDir)
        mkdir(fileDir);
    end
    struct2xml(s,fileId);

end
