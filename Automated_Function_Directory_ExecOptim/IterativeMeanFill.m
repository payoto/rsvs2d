function [varargout]=IterativeMeanFill(varargin)
global IterativeMeanFill_Handle
nOut=nargout(IterativeMeanFill_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=IterativeMeanFill_Handle(varargin{:});
end
