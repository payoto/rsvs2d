function [varargout]=MatchLoops(varargin)
global MatchLoops_Handle
nOut=nargout(MatchLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MatchLoops_Handle(varargin{:});
end
