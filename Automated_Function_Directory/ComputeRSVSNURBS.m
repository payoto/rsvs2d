function [varargout]=ComputeRSVSNURBS(varargin)
% include_NURBSEngine
global ComputeRSVSNURBS_Handle
try
nOut=nargout(ComputeRSVSNURBS_Handle);
catch
include_NURBSEngine
nOut=nargout(ComputeRSVSNURBS_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeRSVSNURBS_Handle(varargin{:});
end
