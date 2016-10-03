function [varargout]=RewriteHistory(varargin)
global RewriteHistory_Handle
nOut=nargout(RewriteHistory_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RewriteHistory_Handle(varargin{:});
end
