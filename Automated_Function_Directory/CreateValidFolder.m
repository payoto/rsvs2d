function [varargout]=CreateValidFolder(varargin)
global CreateValidFolder_Handle
nOut=nargout(CreateValidFolder_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CreateValidFolder_Handle(varargin{:});
end
