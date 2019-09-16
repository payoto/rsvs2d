function [varargout]=ComputeCircularNURBS(varargin)
% include_NURBSEngine
global ComputeCircularNURBS_Handle
nOut=nargout(ComputeCircularNURBS_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeCircularNURBS_Handle(varargin{:});
end
