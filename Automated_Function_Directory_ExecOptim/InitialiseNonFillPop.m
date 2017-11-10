function [varargout]=InitialiseNonFillPop(varargin)
global InitialiseNonFillPop_Handle
nOut=nargout(InitialiseNonFillPop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseNonFillPop_Handle(varargin{:});
end
