function [varargout]=OptimisationOutput_profilepost(varargin)
% OptimisationOutput
global OptimisationOutput_profilepost_Handle
nOut=nargout(OptimisationOutput_profilepost_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationOutput_profilepost_Handle(varargin{:});
end
