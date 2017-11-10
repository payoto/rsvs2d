function [varargout]=GenerateNewPop(varargin)
global GenerateNewPop_Handle
nOut=nargout(GenerateNewPop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateNewPop_Handle(varargin{:});
end
