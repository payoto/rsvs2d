function [varargout]=CutCellFlow(varargin)
global CutCellFlow_Handle
nOut=nargout(CutCellFlow_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=CutCellFlow_Handle(varargin{:});
end
