function [varargout]=EdgeFillInformation(varargin)
global EdgeFillInformation_Handle
nOut=nargout(EdgeFillInformation_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgeFillInformation_Handle(varargin{:});
end
