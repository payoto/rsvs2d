function [varargout]=CompareLoops(varargin)
global CompareLoops_Handle
nOut=nargout(CompareLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompareLoops_Handle(varargin{:});
end
