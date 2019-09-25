function [varargout]=MatchVoltoShape(varargin)
% include_Optimisation
global MatchVoltoShape_Handle
try
nOut=nargout(MatchVoltoShape_Handle);
catch
include_Optimisation
nOut=nargout(MatchVoltoShape_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MatchVoltoShape_Handle(varargin{:});
end
