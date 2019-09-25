function [varargout]=EdgeFillInformationReshape(varargin)
% include_EdgeInformation
global EdgeFillInformationReshape_Handle
try
nOut=nargout(EdgeFillInformationReshape_Handle);
catch
include_EdgeInformation
nOut=nargout(EdgeFillInformationReshape_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgeFillInformationReshape_Handle(varargin{:});
end
