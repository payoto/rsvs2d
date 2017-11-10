function [varargout]=IterativeMinFill(varargin)
global IterativeMinFill_Handle
nOut=nargout(IterativeMinFill_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IterativeMinFill_Handle(varargin{:});
end
