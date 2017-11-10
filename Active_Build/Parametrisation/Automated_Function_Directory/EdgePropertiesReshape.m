function [varargout]=EdgePropertiesReshape(varargin)
global EdgePropertiesReshape_Handle
nOut=nargout(EdgePropertiesReshape_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=EdgePropertiesReshape_Handle(varargin{:});
end
