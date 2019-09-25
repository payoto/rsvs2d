function [varargout]=GenerateEdgeLoop(varargin)
% include_EdgeInformation
global GenerateEdgeLoop_Handle
try
nOut=nargout(GenerateEdgeLoop_Handle);
catch
include_EdgeInformation
nOut=nargout(GenerateEdgeLoop_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateEdgeLoop_Handle(varargin{:});
end
