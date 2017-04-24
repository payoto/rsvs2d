function [varargout]=HandleNonFillVar(varargin)
global HandleNonFillVar_Handle
nOut=nargout(HandleNonFillVar_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=HandleNonFillVar_Handle(varargin{:});
end
