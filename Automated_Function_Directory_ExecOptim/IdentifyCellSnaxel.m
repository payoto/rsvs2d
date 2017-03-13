function [varargout]=IdentifyCellSnaxel(varargin)
global IdentifyCellSnaxel_Handle
nOut=nargout(IdentifyCellSnaxel_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IdentifyCellSnaxel_Handle(varargin{:});
end
