function [varargout]=ComputeErrorVals(varargin)
% include_NURBSEngine
global ComputeErrorVals_Handle
nOut=nargout(ComputeErrorVals_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ComputeErrorVals_Handle(varargin{:});
end
