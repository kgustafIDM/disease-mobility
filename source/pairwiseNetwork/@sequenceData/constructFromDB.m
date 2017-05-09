function obj = constructFromDB(afp,Lineages,strain,strainIdx,varargin)
% constructs sequenceData object from Lineages database
% default behavior is to take all data from a given lineage and strain: ie
%    Lineages.WPV1{1} or Lineages.VDPV2{10} etc
%
% INPUT:
%   afp :: afp database struct or filename. Assumes dates are datenums
%   Lineages :: lineage database struct or filename.  Assumes dates are datenums
%   strain :: strain string ('WPV1', 'WPV3', or 'VDPV2'.  in future, possibly 'Sabin1' etc)
%   strainIdx :: index of substrain in database. only relevant for VDPV2
%                  with multiple emergences. 
%
%   if no imputs, it prompts you to find afp & Lineages Databases, and asks for
%   strain and strain index
%
%   VARARGIN:
%       Default behavior (no additional arguments) is to construct the
%       sequenceData object from all sequences in the lineage.  
%
%       There are multiple options to get subsets of the data.  Options are
%       entered as key:value pairs. If multiple option sets are specified,
%       output is the intersection.
%
%       adminLevel :: gets all sequences in lineage in specified locations
%           examples:
%               'AdminL1',{'loc1','loc2'...}  :: get by standardized location names
%               'AdminL2',[shapeFileId1, shapeFileId2, ...] :: get by shapefile id
%                           shapefilde Ids are currently the only way to
%                           reference location that are guaranteed unambiguous
%
%       dateRange :: get by date range
%           'dateRange', [start1, end1]  possible [start1 end1 start2 end2 ...]
%           'dateUnits', 'days'  or 'years'  (default is 'days')
%       
%       EPID :: get by cellarray of EPIDs
%          'EPID',{'epid1','epid2',...}
%
%       FastaName :: get by cellarray of FastaNames
%           'FastaName', {'fasta1','fasta2',...}
%       
%       Random :: randomly get sequences 
%           'random',[NumRandSeq]  : NumRandSeq = number of sequences wanted
%
%       RandomNormalized :: randomly get sequences with normalized incidence
%           'randomnormalized',[NumRandSeq]  : NumRandSeq = number of sequences wanted
%
%       segment :: get only bases between 2 inclusive limits
%           'segment', [start base position, end position]
%
%       OPTIONAL Provenance data
%           sourceFile :: cell or vector of strings giving sourceFile of lineages DB for provenance
%               'sourceFile' : 'path\file'
%


obj=sequenceData;
options=struct('dateRange',false,'dateUnits','days','EPID',false,'FastaName',false,'AdminL0',false,'AdminL1',false,'AdminL2',false,'ClusterID',false,'random',false,'sourceFile',false,'segment',false,'randomnormalized',false);

% take user input if incomplete arguments given
    if nargin<4
        if nargin<1
            [fileId, filePath] = uigetfile('*afp*.mat','load afp database from PolioInputAnalysis');
            load([filePath,fileId])
            obj.sourceFile=[filePath,fileId];
        end
        if nargin<2
            [fileId, filePath] = uigetfile('*ineage*.mat','load Lineages database from PolioInputAnalysis');
            load([filePath,fileId])
            obj.sourceFile=[obj.sourceFile,'; ',[filePath,fileId]];
        end
        
        if nargin<3,
            strain = input('strain type: ','s');
        end
        
        strainIdx = input('sub-strain index (default = 1): ');
        if isempty(strainIdx)
            strainIdx=1;
        end
        
        moreOpts=input('specify additional options? (y/n) :','s');
        if isempty(moreOpts), moreOpts='n'; end
        if strcmp('y',moreOpts)
            inputs={};
            for F=reshape(fieldnames(options),1,[])
                if ~strcmp(F{:},'sourceFile')
                    tmpInput=input([F{:},': ']);
                    if ~isempty(tmpInput) 
                        inputs{end+1}=F{:};
                        inputs{end+1}=tmpInput;
                    end
                end
            end
            options = keyValuePairVararginHandler(options,inputs);
        end
    end
    strain = regexprep(strain,'''','');
    
% varargin handling
    if nargin>4
        inputs=varargin;
        options = keyValuePairVararginHandler(options,inputs);
    end
    if ~islogical(options.sourceFile)
        obj.soureFile=options.sourceFile;
    end

% store constructor query
    obj.constructorQuery = {strain,num2str(strainIdx)};
    for k=1:length(varargin)
        obj.constructorQuery{end+1}=varargin{k};
    end
    
% get databases if given file names and not DB structs
    switch class(afp)
        case 'char'
            load(afp);    
    end
    switch class(Lineages)
        case 'char'
            load(Lineages);    
    end
    
% get indices into the Lineages and afp databases for the requested sequences
% always starts with all sequneces and then reduces based on options
    LineagesSeqIdx = 1:length(Lineages.(strain){strainIdx}.seq.NT);
    
    [ix, afpSeqIdx] = ismember(Lineages.(strain){strainIdx}.CaseID.EPID(LineagesSeqIdx),afp.CaseID.EPID);
    if any(~ix)
        afpLineagesErrorIdx=find(~ix);
        putvar(afpLineagesErrorIdx);
        warning('sequenceData:afpLineagesMismatch',['EPID mismatch between Lineages and afp databases. Check entries [',num2str(afpLineagesErrorIdx),'] in Lineages.(',strain,'){',num2str(strainIdx),'}.CaseID.EPID .'])
        afpSeqIdx=afpSeqIdx(ix);
    end
    
    % dateRange
        if ~islogical(options.dateRange) && ~isempty(options.dateRange)
            dateIdx=zeros(size(LineagesSeqIdx));
            if mod(length(options.dateRange),2)==0 && any(size(options.dateRange)==1) && any(diff(options.dateRange)>0)
                if strcmp('years',options.dateUnits)
                    options.dateRange = datenumCleaner(options.dateRange,'years');
                else
                    error('dateUnits must either be ''days'' or ''years'' (''days'' by default)')
                end                
                for datePair=reshape(options.dateRange,2,[])
                    dateIdx = dateIdx + (Lineages.(strain){strainIdx}.CaseID.IsolateDate(LineagesSeqIdx)>=datePair(1) & Lineages.(strain){strainIdx}.CaseID.IsolateDate(LineagesSeqIdx)<=datePair(2));
                end
                LineagesSeqIdx = LineagesSeqIdx(logical(dateIdx));
                afpSeqIdx = afpSeqIdx(logical(dateIdx));
            else
                error('sequenceData:dateRangeOption','dateRange option takes pairs in the format [begin1 end1 begin2 end2 ...] only')
            end
        end
    
    % EPID
        if ~islogical(options.EPID) && ~isempty(options.EPID)
            [ix, EPIDIdx]=ismember(options.EPID,Lineages.(strain){strainIdx}.CaseID.EPID(LineagesSeqIdx));
            if any(~ix)
                EPIDErrorIdx=find(~ix);
                putvar(EPIDErrorIdx);
                warning('sequenceData:optionsLineagesMismatch',['Ignorning missing EPIDs. EPID mismatch between options.EPID and Lineages database. Check entries [',num2str(EPIDErrorIdx),'] in options.EPID against Lineages.(',strain,'){',num2str(strainIdx),'}.CaseID.EPID .'])
            end
            LineagesSeqIdx=LineagesSeqIdx(EPIDIdx(ix));
            afpSeqIdx=afpSeqIdx(EPIDIdx(ix));
        end
    
    % FastaName
        if ~islogical(options.FastaName) && ~isempty(options.FastaName)
            [ix, FastaNameIdx]=ismember(options.FastaName,Lineages.(strain){strainIdx}.CaseID.FastaName(LineagesSeqIdx));
            if any(~ix)
                FastaNameErrorIdx=find(~ix);
                putvar(FastaNameErrorIdx);
                warning('sequenceData:optionsLineagesMismatch',['Ignoring missing FastaNames. FastaName mismatch between options.FastaName and Lineages database. Check entries [',num2str(FastaNameErrorIdx),'] in options.FastaName against Lineages.(',strain,'){',num2str(strainIdx),'}.CaseID.FastaName .'])
            end
            LineagesSeqIdx=LineagesSeqIdx(FastaNameIdx(ix));
            afpSeqIdx=afpSeqIdx(FastaNameIdx(ix));
        end    
    
    % AdminL0
        if ~islogical(options.AdminL0) && ~isempty(options.AdminL0)
            AdminL0Idx=zeros(size(afpSeqIdx));
            ix=true(size(options.AdminL0));
            switch class(options.AdminL0)
                case 'cell'
                    for k=1:length(options.AdminL0),
                        tmpIdx=strcmpi(options.AdminL0{k},afp.Location.AdminL0(afpSeqIdx));
                        if all(~tmpIdx), ix(k)=false; end
                        AdminL0Idx = AdminL0Idx + tmpIdx;
                    end
                    AdminL0Idx=logical(AdminL0Idx);
                case 'char'
                    AdminL0Idx = strcmpi(options.AdminL0,afp.Location.AdminL0(afpSeqIdx));
                case 'double'
                    for k=1:length(options.AdminL0)
                        tmpIdx=(options.AdminL0==afp.Location.ShapeFileID.AdminL0(afpSeqIdx));
                        if all(~tmpIdx), ix(k)=false; end
                        AdminL0Idx = AdminL0Idx + tmpIdx;
                    end
                    AdminL0Idx=logical(AdminL0Idx);
            end
            if any(~ix)
                AdminL0ErrorIdx=find(~ix);
                putvar(AdminL0ErrorIdx);
                warning('sequenceData:optionsLineagesMismatch',['Ignoring missing AdminL0. AdminL0 mismatch between options.FastaName and Lineages database. Check entries [',num2str(AdminL0ErrorIdx),'] in options.AdminL0 against afp.Location .'])
            end
            LineagesSeqIdx=LineagesSeqIdx(AdminL0Idx);
            afpSeqIdx=afpSeqIdx(AdminL0Idx);
        end
        
    % AdminL1
        if ~islogical(options.AdminL1) && ~isempty(options.AdminL1)
            AdminL1Idx=zeros(size(afpSeqIdx));
            ix=true(size(options.AdminL1));
            switch class(options.AdminL1)
                case 'cell'
                    for k=1:length(options.AdminL1),
                        tmpIdx=strcmpi(options.AdminL1{k},afp.Location.AdminL1(afpSeqIdx));
                        if all(~tmpIdx), ix(k)=false; end
                        AdminL1Idx = AdminL1Idx + tmpIdx;
                    end
                    AdminL1Idx=logical(AdminL1Idx);
                case 'char'
                    AdminL1Idx = strcmpi(options.AdminL1,afp.Location.AdminL1(afpSeqIdx));
                case 'double'
                    for k=1:length(options.AdminL1)
                        tmpIdx=(options.AdminL1==afp.Location.ShapeFileID.AdminL1(afpSeqIdx));
                        if all(~tmpIdx), ix(k)=false; end
                        AdminL1Idx = AdminL1Idx + tmpIdx;
                    end
                    AdminL1Idx=logical(AdminL1Idx);
            end
            if any(~ix)
                AdminL1ErrorIdx=find(~ix);
                putvar(AdminL1ErrorIdx);
                warning('sequenceData:optionsLineagesMismatch',['Ignoring missing AdminL1. AdminL1 mismatch between options.FastaName and Lineages database. Check entries [',num2str(AdminL1ErrorIdx),'] in options.AdminL1 against afp.Location .'])
            end
            LineagesSeqIdx=LineagesSeqIdx(AdminL1Idx);
            afpSeqIdx=afpSeqIdx(AdminL1Idx);
        end
        
    % AdminL2
        if ~islogical(options.AdminL2) && ~isempty(options.AdminL2)
            AdminL2Idx=zeros(size(afpSeqIdx));
            ix=true(size(options.AdminL2));
            switch class(options.AdminL2)
                case 'cell'
                    for k=1:length(options.AdminL2),
                        tmpIdx=strcmpi(options.AdminL2{k},afp.Location.AdminL2(afpSeqIdx));
                        if all(~tmpIdx), ix(k)=false; end
                        AdminL2Idx = AdminL2Idx + tmpIdx;
                    end
                    AdminL2Idx=logical(AdminL2Idx);
                case 'char'
                    AdminL2Idx = strcmpi(options.AdminL2,afp.Location.AdminL2(afpSeqIdx));
                case 'double'
                    for k=1:length(options.AdminL2)
                        tmpIdx=(options.AdminL2==afp.Location.ShapeFileID.AdminL2(afpSeqIdx));
                        if all(~tmpIdx), ix(k)=false; end
                        AdminL2Idx = AdminL2Idx + tmpIdx;
                    end
                    AdminL2Idx=logical(AdminL2Idx);
            end
            if any(~ix)
                AdminL2ErrorIdx=find(~ix);
                putvar(AdminL2ErrorIdx);
                warning('sequenceData:optionsLineagesMismatch',['Ignoring missing AdminL2. AdminL2 mismatch between options.FastaName and Lineages database. Check entries [',num2str(AdminL2ErrorIdx),'] in options.AdminL2 against afp.Location .'])
            end
            LineagesSeqIdx=LineagesSeqIdx(AdminL2Idx);
            afpSeqIdx=afpSeqIdx(AdminL2Idx);
        end

    % ClusterID
        if ~islogical(options.ClusterID) && ~isempty(options.ClusterID)
            ClusterIDIdx=zeros(size(afpSeqIdx));
            ix=true(size(options.ClusterID));
            switch class(options.ClusterID)
                case 'cell'
                    for k=1:length(options.ClusterID),
                        tmpIdx=strcmpi(options.ClusterID{k},afp.Location.ClusterID(afpSeqIdx));
                        if all(~tmpIdx), ix(k)=false; end
                        ClusterIDIdx = ClusterIDIdx + tmpIdx;
                    end
                    ClusterIDIdx=logical(ClusterIDIdx);
                case 'char'
                    ClusterIDIdx = strcmpi(options.ClusterID,afp.Location.ClusterID(afpSeqIdx));
                case 'double'
                    for k=1:length(options.ClusterID)
                        tmpIdx=(options.ClusterID==afp.Location.ShapeFileID.ClusterID(afpSeqIdx));
                        if all(~tmpIdx), ix(k)=false; end
                        ClusterIDIdx = ClusterIDIdx + tmpIdx;
                    end
                    ClusterIDIdx=logical(ClusterIDIdx);
            end
            if any(~ix)
                ClusterIDErrorIdx=find(~ix);
                putvar(ClusterIDErrorIdx);
                warning('sequenceData:optionsLineagesMismatch',['Ignoring missing ClusterID. ClusterID mismatch between options.FastaName and Lineages database. Check entries [',num2str(ClusterIDErrorIdx),'] in options.ClusterID against afp.Location .'])
            end
            LineagesSeqIdx=LineagesSeqIdx(ClusterIDIdx);
            afpSeqIdx=afpSeqIdx(ClusterIDIdx);
        end
        
    % random selection
        if ~islogical(options.random) && ~isempty(options.random) && islogical(options.randomnormalized)
            if length(LineagesSeqIdx)<=options.random
                warning('Ignoring random option: number of available sequences after other options are processed is less than or equal to the number of sequences requested.')
            else
                randIdx=randperm(length(LineagesSeqIdx),options.random);
                LineagesSeqIdx=LineagesSeqIdx(randIdx);
                afpSeqIdx=afpSeqIdx(randIdx);                
            end
        elseif ~islogical(options.randomnormalized) && ~isempty(options.randomnormalized) && islogical(options.random)
            if length(LineagesSeqIdx)<=options.random
                warning('Ignoring random option: number of available sequences after other options are processed is less than or equal to the number of sequences requested.')
            else
                numPerBin=ceil(options.randomnormalized/20);
                y=linspace(min(Lineages.(strain){strainIdx}.CaseID.IsolateDate),max(Lineages.(strain){strainIdx}.CaseID.IsolateDate),21);
                y=y(2:end);
                randIdx=zeros(options.randomnormalized,1);
                for qq=1:ceil(options.randomnormalized/20)
                    if qq==1
                        count=sum(Lineages.(strain){strainIdx}.CaseID.IsolateDate<=y(1));
                        if count>numPerBin
                            randIdx(1:numPerBin)=randperm(count,numPerBin);
                        else
                            randIdx(1:count)=1:count;
                        end
                    else
                        count=sum(Lineages.(strain){strainIdx}.CaseID.IsolateDate>y(qq-1) & Lineages.(strain){strainIdx}.CaseID.IsolateDate<=y(qq));
                        if count>numPerBin
                            randIdx(((qq-1)*numPerBin +1 :qq*numPerBin))=sum(Lineages.(strain){strainIdx}.CaseID.IsolateDate<=y(qq-1))+randperm(count,numPerBin);
                        else
                            randIdx(((qq-1)*numPerBin +1 :(qq-1)*numPerBin+count))=sum(Lineages.(strain){strainIdx}.CaseID.IsolateDate<=y(qq-1))+(1:count);
                        end
                    end
                end
                randIdx=randIdx(randIdx>0);
                LineagesSeqIdx=LineagesSeqIdx(randIdx);
                afpSeqIdx=afpSeqIdx(randIdx);                
            end
        else
            warning('can either use random or randomnormalized. assuming randomnormalized only');
        end
        
        
% check to see if there is any data to get after all the options are processed
    if isempty(LineagesSeqIdx)
        error('sequenceData:noSequencesFound','there are no cases consistent with all of your options.')
    end

% get data
    numSeq=length(LineagesSeqIdx);
    obj.isolateDate = datenumCleaner(Lineages.(strain){strainIdx}.CaseID.IsolateDate(LineagesSeqIdx));       
    obj.EPID = Lineages.(strain){strainIdx}.CaseID.EPID(LineagesSeqIdx)';
    obj.FastaName = Lineages.(strain){strainIdx}.CaseID.FastaName(LineagesSeqIdx)';
    obj.strain = repmat({strain},numSeq,1);
    obj.phylogeneticDomain = repmat({Lineages.(strain){strainIdx}.PhylogeneticDomain},numSeq,1);
    
    if isfield(Lineages.(strain){strainIdx}.CaseID,'ClusterID')
        obj.ClusterID = Lineages.(strain){strainIdx}.CaseID.ClusterID(LineagesSeqIdx)';
    end

    obj.sequence = cell(numSeq,1);
    
    for k=1:numSeq,
        obj.sequence{k} = upper(Lineages.(strain){strainIdx}.seq.NT{LineagesSeqIdx(k)}(Lineages.(strain){strainIdx}.seq.CDSrange(1,LineagesSeqIdx(k)):Lineages.(strain){strainIdx}.seq.CDSrange(2,LineagesSeqIdx(k))));
    end
    % are all sequences the same length?  If not, take global alignment with all its gap ('-') data
        seqLen=cellfun(@length,obj.sequence);
        if any(~(seqLen==seqLen(1)))
            warning('sequences are of variable length')
            for k=1:numSeq,
                obj.sequence{k}=upper(Lineages.(strain){strainIdx}.GlobalAlignment(LineagesSeqIdx(k),:));
            end
        end
    
    % take only a segment
        if ~islogical(options.segment)
            for k=1:numSeq,
                obj.sequence{k}=obj.sequence{k}(options.segment(1):options.segment(2));
            end
        end
    
    obj.AdminL0 = lower(afp.Location.AdminL0(afpSeqIdx)');
    obj.AdminL1 = lower(afp.Location.AdminL1(afpSeqIdx)');
    obj.AdminL2 = lower(afp.Location.AdminL2(afpSeqIdx)');
    if ~all(isnan(afp.Location.ShapeFileID.AdminL0(afpSeqIdx)))  % one country at a time right now, and this field is typically NaN. I would rather pass nothing than pass NaN
        obj.shapeFileIdAdminL0=afp.Location.ShapeFileID.AdminL0(afpSeqIdx)';
    end
    obj.shapeFileIdAdminL1=afp.Location.ShapeFileID.AdminL1(afpSeqIdx)';
    obj.shapeFileIdAdminL2=afp.Location.ShapeFileID.AdminL2(afpSeqIdx)';
    
% clean up    
    checkErrors(obj)
end