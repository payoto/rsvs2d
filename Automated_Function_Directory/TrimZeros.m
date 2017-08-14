function [varargout]=TrimZeros(varargin)
% include_Utilities
global TrimZeros_Handle
nOut=nargout(TrimZeros_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TrimZeros_Handle(varargin{:});
end
