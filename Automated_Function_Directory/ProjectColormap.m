function [varargout]=ProjectColorMap(varargin)
% include_PostProcessing
global ProjectColorMap_Handle
try
nOut=nargout(ProjectColorMap_Handle);
catch
include_PostProcessing
nOut=nargout(ProjectColorMap_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProjectColorMap_Handle(varargin{:});
end
