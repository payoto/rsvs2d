function [varargout]=ComputeRSVSNURBS(varargin)
% include_NURBSEngine
global ComputeRSVSNURBS_Handle
nOut=nargout(ComputeRSVSNURBS_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeRSVSNURBS_Handle(varargin{:});
end
