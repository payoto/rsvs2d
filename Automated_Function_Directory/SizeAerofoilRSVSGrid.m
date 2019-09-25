function [varargout]=SizeAerofoilRSVSGrid(varargin)
% include_Utilities
global SizeAerofoilRSVSGrid_Handle
try
nOut=nargout(SizeAerofoilRSVSGrid_Handle);
catch
include_Utilities
nOut=nargout(SizeAerofoilRSVSGrid_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SizeAerofoilRSVSGrid_Handle(varargin{:});
end
