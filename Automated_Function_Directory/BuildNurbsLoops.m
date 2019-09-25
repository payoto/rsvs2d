function [varargout]=BuildNurbsLoops(varargin)
% include_NURBSEngine
global BuildNurbsLoops_Handle
try
nOut=nargout(BuildNurbsLoops_Handle);
catch
include_NURBSEngine
nOut=nargout(BuildNurbsLoops_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildNurbsLoops_Handle(varargin{:});
end
