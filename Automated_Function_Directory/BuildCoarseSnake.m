function [varargout]=BuildCoarseSnake(varargin)
% include_NURBSEngine
global BuildCoarseSnake_Handle
try
nOut=nargout(BuildCoarseSnake_Handle);
catch
include_NURBSEngine
nOut=nargout(BuildCoarseSnake_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildCoarseSnake_Handle(varargin{:});
end
