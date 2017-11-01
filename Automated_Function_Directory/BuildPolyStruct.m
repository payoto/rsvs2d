function [varargout]=BuildPolyStruct(varargin)
% include_PostProcessing
global BuildPolyStruct_Handle
nOut=nargout(BuildPolyStruct_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildPolyStruct_Handle(varargin{:});
end
