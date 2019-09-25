function [varargout]=FinishLoops(varargin)
% include_SnakeParam
global FinishLoops_Handle
try
nOut=nargout(FinishLoops_Handle);
catch
include_SnakeParam
nOut=nargout(FinishLoops_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FinishLoops_Handle(varargin{:});
end
