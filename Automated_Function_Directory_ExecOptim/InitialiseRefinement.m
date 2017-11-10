function [varargout]=InitialiseRefinement(varargin)
global InitialiseRefinement_Handle
nOut=nargout(InitialiseRefinement_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InitialiseRefinement_Handle(varargin{:});
end
