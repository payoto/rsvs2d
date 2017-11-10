function [varargout]=OptimisationOutput_iteration(varargin)
% OptimisationOutput
global OptimisationOutput_iteration_Handle
nOut=nargout(OptimisationOutput_iteration_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationOutput_iteration_Handle(varargin{:});
end
