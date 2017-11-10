function [varargout]=FindProfile(varargin)
% OptimisationOutput
global FindProfile_Handle
nOut=nargout(FindProfile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindProfile_Handle(varargin{:});
end
