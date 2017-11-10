function [varargout]=OptimisationOutput_iterationFullPop(varargin)
% OptimisationOutput
global OptimisationOutput_iterationFullPop_Handle
nOut=nargout(OptimisationOutput_iterationFullPop_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationOutput_iterationFullPop_Handle(varargin{:});
end
