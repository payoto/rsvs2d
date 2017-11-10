function [varargout]=CCWLoop(varargin)
% include_SnakeParam
global CCWLoop_Handle
nOut=nargout(CCWLoop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CCWLoop_Handle(varargin{:});
end
