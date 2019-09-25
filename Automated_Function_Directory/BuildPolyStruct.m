function [varargout]=BuildPolyStruct(varargin)
% include_PostProcessing
global BuildPolyStruct_Handle
try
nOut=nargout(BuildPolyStruct_Handle);
catch
include_PostProcessing
nOut=nargout(BuildPolyStruct_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildPolyStruct_Handle(varargin{:});
end
