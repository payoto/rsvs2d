function [varargout]=EdgePropertiesReshape(varargin)
% include_EdgeInformation
global EdgePropertiesReshape_Handle
try
nOut=nargout(EdgePropertiesReshape_Handle);
catch
include_EdgeInformation
nOut=nargout(EdgePropertiesReshape_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgePropertiesReshape_Handle(varargin{:});
end
