function [varargout]=CalculateEdgeLengths(varargin)
% include_EdgeInformation
global CalculateEdgeLengths_Handle
try
nOut=nargout(CalculateEdgeLengths_Handle);
catch
include_EdgeInformation
nOut=nargout(CalculateEdgeLengths_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CalculateEdgeLengths_Handle(varargin{:});
end
