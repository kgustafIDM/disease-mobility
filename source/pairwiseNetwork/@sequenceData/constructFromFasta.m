function obj = constructFromFasta(obj,fileId,varargin)
% constructs sequenceData object from Lineages database
% default behavior is to take all data from a given lineage and strain: ie
%    Lineages.WPV1{1} or Lineages.VDPV2{10} etc
%
% INPUT:
%   soureFile : 'fileIdString' pointing to source fasta-formatted file
%
%   VARARGIN:
%       Default behavior (no additional arguments) is to construct the
%       sequenceData object from all sequences in the lineage.  
%
%       There are multiple options to get subsets of the data.  Options are
%       entered as key:value pairs. If multiple option sets are specified,
%       output is the intersection.
%
%       soureFormat :: parses headers according to known source format
%           'sourceFormat' , 'genbank' (default) also 'EMBL','DDBJ' and 'generic'
%       
%       Random :: randomly get sequences 
%           'random',[NumRandSeq]  : NumRandSeq = number of sequences wanted
%
%       segment :: get only bases between 2 inclusive limits
%           'segment', [start base position, end position]
%
% TO DO: ignoreAttributes is currently dangerous, as sequence and trait
% info links are easily broken

options=struct('sourceFormat','genbank','random',false,'segment',false,'ignoreAttributes',false);
obj.sourceFile=fileId;

% varargin handling
    if nargin>1
        inputs=varargin;
        options = keyValuePairVararginHandler(options,inputs);
    end

% store constructor query
    obj.constructorQuery = {};
    for k=1:length(varargin)
        obj.constructorQuery{end+1}=varargin{k};
    end
    
% load file
    fastaIn=fastaread(fileId);
    numSeq=length(fastaIn);
    seqIdx=1:numSeq;

% random selection
    if ~islogical(options.random) && ~isempty(options.random)
        if numSeq<=options.random
            warning('Ignoring random option: number of available sequences after other options are processed is less than or equal to the number of sequences requested.')
        else
            seqIdx=sort(randperm(numSeq,options.random));
        end
        numSeq=length(seqIdx);
    end    
    
    
% store sequences
    obj.sequence=cell(numSeq,1);
    for k=1:numSeq
        obj.sequence{k}=fastaIn(seqIdx(k)).Sequence;
    end
    % are all sequences the same length?  If not, take global alignment with all its gap ('-') data
        seqLen=cellfun(@length,obj.sequence);
        if any(~(seqLen==seqLen(1)))
            warning('sequences are of variable length')
        end
    
    % take only a segment
        if ~islogical(options.segment)
            for k=1:numSeq,
                obj.sequence{k}=obj.sequence{k}(options.segment(1):options.segment(2));
            end
        end    
    
% shoehorn header properties into attribute fields
    if ~options.ignoreAttributes
        obj.FastaName=cell(numSeq,1);
        obj.attributes=cell(numSeq,1);
        switch options.sourceFormat
            case {'genbank','EMBL','DDBJ'}
                obj.EPID=cell(numSeq,1);
                for k=1:numSeq
                    tmp=get_tokens(fastaIn(seqIdx(k)).Header,'|');
                    obj.EPID{k}=[tmp{1},'-',tmp{2}];
                    obj.FastaName{k}=[tmp{3},'-',tmp{4}];
                    for m=5:length(tmp)
                        obj.attributes{k}=[obj.attributes{k},'-',tmp{m}];
                    end
                end
            case 'generic'
                for k=1:numSeq,
                    bin=length(fastaIn(seqIdx(k)).Header);
                    if bin>20,
                        bin=20;
                    end
                    obj.FastaName{k}=fastaIn(seqIdx(k)).Header(1:bin);
                    obj.attributes{k}=fastaIn(seqIdx(k)).Header;
                end
            otherwise
            % request user input    
        end
    end
end