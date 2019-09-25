function [varargout]=ComputeCircularNURBS(varargin)
% include_NURBSEngine
global ComputeCircularNURBS_Handle
try
nOut=nargout(ComputeCircularNURBS_Handle);
catch
include_NURBSEngine
nOut=nargout(ComputeCircularNURBS_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeCircularNURBS_Handle(varargin{:});
end
