function [varargout]=ProjectColormap(varargin)
% include_PostProcessing
global ProjectColormap_Handle
try
nOut=nargout(ProjectColormap_Handle);
catch
include_PostProcessing
nOut=nargout(ProjectColormap_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProjectColormap_Handle(varargin{:});
end
