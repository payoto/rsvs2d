function [varargout]=MatchVoltoShape(varargin)
global MatchVoltoShape_Handle
nOut=nargout(MatchVoltoShape_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MatchVoltoShape_Handle(varargin{:});
end
