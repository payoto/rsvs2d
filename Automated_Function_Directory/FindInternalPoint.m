function [varargout]=FindInternalPoint(varargin)
% include_PostProcessing
global FindInternalPoint_Handle
nOut=nargout(FindInternalPoint_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindInternalPoint_Handle(varargin{:});
end
