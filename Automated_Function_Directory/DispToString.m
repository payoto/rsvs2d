function [varargout]=DispToString(varargin)
% include_PostProcessing
global DispToString_Handle
try
nOut=nargout(DispToString_Handle);
catch
include_PostProcessing
nOut=nargout(DispToString_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DispToString_Handle(varargin{:});
end
