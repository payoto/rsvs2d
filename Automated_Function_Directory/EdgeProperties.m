function [varargout]=EdgeProperties(varargin)
% include_EdgeInformation
global EdgeProperties_Handle
nOut=nargout(EdgeProperties_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgeProperties_Handle(varargin{:});
end
