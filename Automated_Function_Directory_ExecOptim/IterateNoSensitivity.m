function [varargout]=IterateNoSensitivity(varargin)
global IterateNoSensitivity_Handle
nOut=nargout(IterateNoSensitivity_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IterateNoSensitivity_Handle(varargin{:});
end
