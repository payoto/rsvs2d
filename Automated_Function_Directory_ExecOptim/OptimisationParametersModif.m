function [varargout]=OptimisationParametersModif(varargin)
global OptimisationParametersModif_Handle
nOut=nargout(OptimisationParametersModif_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationParametersModif_Handle(varargin{:});
end
