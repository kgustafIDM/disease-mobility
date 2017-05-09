function obj = combine(varargin)

if ~all(strcmp('sequenceData',cellfun(@class,varargin,'uniformOutput',false)))
    error('all inputs must be sequenceData objects')
end

numSeq=0;
for k=1:length(varargin)
    numSeq=numSeq+varargin{k}.numSeq;
end

obj=sequenceData;

%pre-declare
    obj.sourceFile =[];
    obj.constructorQuery =[];
    obj.attributes ={};

count=0;
for k=1:length(varargin)
    obj.sourceFile=[obj.sourceFile,'; ',varargin{k}.sourceFile];
    obj.constructorQuery=[obj.constructorQuery,'; ',varargin{k}.constructorQuery];
    obj.attributes{end+1}=varargin{k}.attributes;
    for m=1:varargin{k}.numSeq,
        count=count+1;
        for name = reshape(fieldnames(varargin{k}),1,[])
            inpName=name{1};
            if any(strcmp(inpName,obj.cellFields)) 
                if count==1 ,
                    obj.(inpName)=cell(numSeq,1);
                end
                if ~isempty(varargin{k}.(inpName))
                    obj.(inpName){count}=varargin{k}.(inpName){m};
                end
                if count==numSeq,
                    if all(cellfun(@isempty,obj.(inpName))),
                        obj.(inpName)={};
                    end
                end
            elseif any(strcmp(inpName,obj.doubleFields)),
                if count==1,
                    obj.(inpName)=zeros(numSeq,1);
                end
                if ~isempty(varargin{k}.(inpName))
                    obj.(inpName)(count)=varargin{k}.(inpName)(m);
                end
                if count==numSeq,
                    if all(obj.(inpName)==0),
                        obj.(inpName)=[];
                    end
                end
            end
        end
                
    end
end
checkErrors(obj)
end