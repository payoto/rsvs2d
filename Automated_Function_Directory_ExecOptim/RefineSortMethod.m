function [varargout]=RefineSortMethod(varargin)
global RefineSortMethod_Handle
nOut=nargout(RefineSortMethod_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RefineSortMethod_Handle(varargin{:});
end
