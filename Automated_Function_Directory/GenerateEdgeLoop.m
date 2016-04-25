function [varargout]=GenerateEdgeLoop(varargin)
global GenerateEdgeLoop_Handle
nOut=nargout(GenerateEdgeLoop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateEdgeLoop_Handle(varargin{:});
end
