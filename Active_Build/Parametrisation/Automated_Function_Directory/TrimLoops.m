function [varargout]=TrimLoops(varargin)
global TrimLoops_Handle
nOut=nargout(TrimLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TrimLoops_Handle(varargin{:});
end
