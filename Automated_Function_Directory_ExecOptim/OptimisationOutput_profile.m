function [varargout]=OptimisationOutput_profile(varargin)
% OptimisationOutput
global OptimisationOutput_profile_Handle
nOut=nargout(OptimisationOutput_profile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimisationOutput_profile_Handle(varargin{:});
end
