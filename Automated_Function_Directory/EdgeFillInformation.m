function [varargout]=EdgeFillInformation(varargin)
% include_EdgeInformation
global EdgeFillInformation_Handle
try
nOut=nargout(EdgeFillInformation_Handle);
catch
include_EdgeInformation
nOut=nargout(EdgeFillInformation_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgeFillInformation_Handle(varargin{:});
end
