function [varargout]=TrimZerosImage(varargin)
% include_Utilities
global TrimZerosImage_Handle
nOut=nargout(TrimZerosImage_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TrimZerosImage_Handle(varargin{:});
end
