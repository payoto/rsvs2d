function [varargout]=HandleRefinement(varargin)
global HandleRefinement_Handle
nOut=nargout(HandleRefinement_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=HandleRefinement_Handle(varargin{:});
end
