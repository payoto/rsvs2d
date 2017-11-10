function [varargout]=BuildDefaultFillRef(varargin)
global BuildDefaultFillRef_Handle
nOut=nargout(BuildDefaultFillRef_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildDefaultFillRef_Handle(varargin{:});
end
