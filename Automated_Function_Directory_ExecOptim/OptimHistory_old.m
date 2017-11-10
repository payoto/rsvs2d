function [varargout]=OptimHistory_old(varargin)
% OptimisationOutput
global OptimHistory_old_Handle
nOut=nargout(OptimHistory_old_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimHistory_old_Handle(varargin{:});
end
