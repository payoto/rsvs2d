function [varargout]=FELINESEGHeader(varargin)
global FELINESEGHeader_Handle
nOut=nargout(FELINESEGHeader_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FELINESEGHeader_Handle(varargin{:});
end
