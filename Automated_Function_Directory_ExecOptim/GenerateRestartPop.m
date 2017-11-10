function [varargout]=GenerateRestartPop(varargin)
global GenerateRestartPop_Handle
nOut=nargout(GenerateRestartPop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateRestartPop_Handle(varargin{:});
end
