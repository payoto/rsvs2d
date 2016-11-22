function [a]=opstruct(operator,varargin)
    % Takes in an operator and then returns
    
    if numel(varargin)==2
        a=varargin{1};
        b=varargin{2};
        operator2='';
        if ~isstruct(a) && ~isstruct(b)
            error('Function not defined for nonstructure arguments')
        elseif ~isstruct(a)
            temp=a;a=b;b=temp;
            if strcmp(operator,'-');operator2='-';end
            if ~isempty(regexp(operator,'/', 'once'));operator2=['1',operator];end
            if ~isempty(regexp(operator,'<', 'once'));operator=regexprep('<','<','>');end
            if ~isempty(regexp(operator,'>', 'once'));operator=regexprep('>','>','<');end
        end

        fA=fieldnames(a);
        if isstruct(b)
            if isStructEqual(a,b)
                for ii=1:numel(fA)
                    a.(fA{ii})=eval([operator2,'(a.(fA{ii})',operator,'b.(fA{ii}))']);
                end
            else
                error('Unequal structures a and b')
            end
        else
            for ii=1:numel(fA)
                a.(fA{ii})=eval([operator2,'(a.(fA{ii})',operator,'b)']);
            end
        end
    else
        isStructArgIn=cellfun(@isstruct,varargin);
        structArgIn=find(isStructArgIn);
        if isempty(structArgIn)
            error('No structures in the input')
        end
        for ii=2:numel(structArgIn)
            if ~isStructEqual(varargin{structArgIn(1)},varargin{structArgIn(ii)})
                error('Unequal structures %i and 1',ii)
            end
        end
        a=varargin{structArgIn(1)}(1);
        fA=fieldnames(a);
        opstring='(';
        for ii=1:numel(varargin)
            opstring=[opstring,'varargin{',int2str(ii),'}'];
            if isStructArgIn(ii)
                opstring=[opstring,'.(fA{ii}),'];
            else
                opstring=[opstring,','];
            end
        end
        opstring(end)=')';
        
        if ischar(operator)
            for ii=1:numel(fA)
                a.(fA{ii})=eval([operator,opstring]);
            end
        else
            for ii=1:numel(fA)
                a.(fA{ii})=eval(['operator',opstring]);
            end
        end
        
    end
end


function [a]=opstruct2(a,b,operator)
    
    fA=fieldnames(a);
    operator2='';
    if ~isstruct(a) && ~isstruct(b)
        error('Function not defined for nonstructure arguments')
    elseif ~isstruct(a)
        temp=a;a=b;b=temp;
        if strcmp(operator,'-');operator2='-';end
        if ~isempty(regexp(operator,'/', 'once'));operator2=['1',operator];end
    end
        
    if isstruct(b)
        if isStructEqual(a,b)
            for ii=1:numel(fA)
                a.(fA{ii})=eval([operator2,'(a.(fA{ii})',operator,'b.(fA{ii}))']);
            end
        else
            error('Unequal structures a and b')
        end
    else
        for ii=1:numel(fA)
            a.(fA{ii})=eval([operator2,'(a.(fA{ii})',operator,'b)']);
        end
    end
end

