function [varargout]=SizeAerofoilRSVSGrid(varargin)
% include_Utilities
global SizeAerofoilRSVSGrid_Handle
nOut=nargout(SizeAerofoilRSVSGrid_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=SizeAerofoilRSVSGrid_Handle(varargin{:});
end
