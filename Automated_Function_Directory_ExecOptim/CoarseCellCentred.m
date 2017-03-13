function [varargout]=CoarseCellCentred(varargin)
global CoarseCellCentred_Handle
nOut=nargout(CoarseCellCentred_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CoarseCellCentred_Handle(varargin{:});
end
