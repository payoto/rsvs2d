function [varargout]=ComputeCircularNURBS2(varargin)
% include_NURBSEngine
global ComputeCircularNURBS2_Handle
try
nOut=nargout(ComputeCircularNURBS2_Handle);
catch
include_NURBSEngine
nOut=nargout(ComputeCircularNURBS2_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeCircularNURBS2_Handle(varargin{:});
end
