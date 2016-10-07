function [varargout]=AddZeroLayer(varargin)
global AddZeroLayer_Handle
nOut=nargout(AddZeroLayer_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=AddZeroLayer_Handle(varargin{:});
end
