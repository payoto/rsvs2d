function [varargout]=IdentifyUniqueOptions(varargin)
% include_PostProcessing
global IdentifyUniqueOptions_Handle
try
nOut=nargout(IdentifyUniqueOptions_Handle);
catch
include_PostProcessing
nOut=nargout(IdentifyUniqueOptions_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IdentifyUniqueOptions_Handle(varargin{:});
end
