function [varargout]=BuildMatrix(varargin)
global BuildMatrix_Handle
nOut=nargout(BuildMatrix_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildMatrix_Handle(varargin{:});
end
