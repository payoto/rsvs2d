function [varargout]=OptimHistory(varargin)
% OptimisationOutput
global OptimHistory_Handle
nOut=nargout(OptimHistory_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OptimHistory_Handle(varargin{:});
end
