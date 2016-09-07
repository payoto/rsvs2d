function [varargout]=CalculatePolyArea(varargin)
global CalculatePolyArea_Handle
nOut=nargout(CalculatePolyArea_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculatePolyArea_Handle(varargin{:});
end
