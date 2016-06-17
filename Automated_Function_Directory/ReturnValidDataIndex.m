function [varargout]=ReturnValidDataIndex(varargin)
global ReturnValidDataIndex_Handle
nOut=nargout(ReturnValidDataIndex_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReturnValidDataIndex_Handle(varargin{:});
end
