function [varargout]=EdgeFillInformationReshape(varargin)
% include_EdgeInformation
global EdgeFillInformationReshape_Handle
nOut=nargout(EdgeFillInformationReshape_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgeFillInformationReshape_Handle(varargin{:});
end
