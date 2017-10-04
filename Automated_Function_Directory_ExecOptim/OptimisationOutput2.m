function [varargout]=OptimisationOutput2(varargin)
% OptimisationOutput
global OptimisationOutput2_Handle
nOut=nargout(OptimisationOutput2_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationOutput2_Handle(varargin{:});
end
