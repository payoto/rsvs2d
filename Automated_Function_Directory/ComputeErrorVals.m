function [varargout]=ComputeErrorVals(varargin)
% include_NURBSEngine
global ComputeErrorVals_Handle
try
nOut=nargout(ComputeErrorVals_Handle);
catch
include_NURBSEngine
nOut=nargout(ComputeErrorVals_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeErrorVals_Handle(varargin{:});
end
