function [varargout]=BuildCoarseSnake(varargin)
% include_NURBSEngine
global BuildCoarseSnake_Handle
nOut=nargout(BuildCoarseSnake_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildCoarseSnake_Handle(varargin{:});
end
