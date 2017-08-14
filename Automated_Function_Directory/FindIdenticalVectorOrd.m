function [varargout]=FindIdenticalVectorOrd(varargin)
% include_Utilities
global FindIdenticalVectorOrd_Handle
nOut=nargout(FindIdenticalVectorOrd_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindIdenticalVectorOrd_Handle(varargin{:});
end
