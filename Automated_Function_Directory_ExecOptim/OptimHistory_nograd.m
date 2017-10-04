function [varargout]=OptimHistory_nograd(varargin)
% OptimisationOutput
global OptimHistory_nograd_Handle
nOut=nargout(OptimHistory_nograd_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimHistory_nograd_Handle(varargin{:});
end
