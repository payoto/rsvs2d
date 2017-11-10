function [varargout]=OptimHistory_grad(varargin)
% OptimisationOutput
global OptimHistory_grad_Handle
nOut=nargout(OptimHistory_grad_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimHistory_grad_Handle(varargin{:});
end
