function [varargout]=LoopToGridFill(varargin)
% include_Optimisation
global LoopToGridFill_Handle
try
nOut=nargout(LoopToGridFill_Handle);
catch
include_Optimisation
nOut=nargout(LoopToGridFill_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LoopToGridFill_Handle(varargin{:});
end
