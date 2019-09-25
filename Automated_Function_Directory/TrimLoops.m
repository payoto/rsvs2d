function [varargout]=TrimLoops(varargin)
% include_PostProcessing
global TrimLoops_Handle
try
nOut=nargout(TrimLoops_Handle);
catch
include_PostProcessing
nOut=nargout(TrimLoops_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TrimLoops_Handle(varargin{:});
end
