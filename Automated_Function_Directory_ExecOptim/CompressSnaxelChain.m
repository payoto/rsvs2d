function [varargout]=CompressSnaxelChain(varargin)
global CompressSnaxelChain_Handle
nOut=nargout(CompressSnaxelChain_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CompressSnaxelChain_Handle(varargin{:});
end
