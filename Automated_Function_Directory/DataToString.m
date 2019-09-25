function [varargout]=DataToString(varargin)
% include_PostProcessing
global DataToString_Handle
try
nOut=nargout(DataToString_Handle);
catch
include_PostProcessing
nOut=nargout(DataToString_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DataToString_Handle(varargin{:});
end
