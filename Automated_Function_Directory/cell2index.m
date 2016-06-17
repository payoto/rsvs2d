function [varargout]=cell2index(varargin)
global cell2index_Handle
nOut=nargout(cell2index_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=cell2index_Handle(varargin{:});
end
