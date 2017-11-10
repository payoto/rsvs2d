function [varargout]=FinishLoops(varargin)
% include_SnakeParam
global FinishLoops_Handle
nOut=nargout(FinishLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FinishLoops_Handle(varargin{:});
end
