function [varargout]=FindKnownOptimCutCell(varargin)
global FindKnownOptimCutCell_Handle
nOut=nargout(FindKnownOptimCutCell_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindKnownOptimCutCell_Handle(varargin{:});
end
