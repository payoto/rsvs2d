function [varargout]=CalculateEdgeLengths(varargin)
global CalculateEdgeLengths_Handle
nOut=nargout(CalculateEdgeLengths_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculateEdgeLengths_Handle(varargin{:});
end
