function [varargout]=MakePathCompliant(varargin)
% include_PostProcessing
global MakePathCompliant_Handle
try
nOut=nargout(MakePathCompliant_Handle);
catch
include_PostProcessing
nOut=nargout(MakePathCompliant_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakePathCompliant_Handle(varargin{:});
end
