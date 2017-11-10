function [varargout]=AnisotropicRefinement(varargin)
global AnisotropicRefinement_Handle
nOut=nargout(AnisotropicRefinement_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=AnisotropicRefinement_Handle(varargin{:});
end
