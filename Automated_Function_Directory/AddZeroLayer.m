function [varargout]=AddZeroLayer(varargin)
% include_Utilities
global AddZeroLayer_Handle
try
nOut=nargout(AddZeroLayer_Handle);
catch
include_Utilities
nOut=nargout(AddZeroLayer_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=AddZeroLayer_Handle(varargin{:});
end
