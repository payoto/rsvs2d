function [varargout]=ReFracGrids(varargin)
global ReFracGrids_Handle
nOut=nargout(ReFracGrids_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ReFracGrids_Handle(varargin{:});
end
