function [varargout]=BuildIterRes(varargin)
% OptimisationOutput
global BuildIterRes_Handle
nOut=nargout(BuildIterRes_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=BuildIterRes_Handle(varargin{:});
end
