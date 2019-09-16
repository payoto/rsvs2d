function [varargout]=LoopToGridFill(varargin)
% include_Optimisation
global LoopToGridFill_Handle
nOut=nargout(LoopToGridFill_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LoopToGridFill_Handle(varargin{:});
end
