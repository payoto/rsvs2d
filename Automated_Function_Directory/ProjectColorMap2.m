function [varargout]=ProjectColorMap2(varargin)
% include_PostProcessing
global ProjectColorMap2_Handle
try
nOut=nargout(ProjectColorMap2_Handle);
catch
include_PostProcessing
nOut=nargout(ProjectColorMap2_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProjectColorMap2_Handle(varargin{:});
end
