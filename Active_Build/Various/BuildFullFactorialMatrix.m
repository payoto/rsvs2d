function [ffmat]=BuildFullFactorialMatrix(varargin)
   % Builds a full factorial matrix for any number of numeric arrays
    
    ffSize=[prod(cellfun(@numel,varargin)),nargin];
    
    ffmat=zeros(ffSize);
    
    freq=1;
    for ii=1:nargin
        temp=reshape(varargin{ii},[numel(varargin{ii}),1]);
        ffmat(:,ii)=temp(mod(ceil((1:ffSize(1))/freq)-1,numel(varargin{ii}))+1);
        freq=freq*numel(varargin{ii});
    end
    
    
end