function [varargout]=MatchNURBSpatchArea(varargin)
% include_NURBSEngine
global MatchNURBSpatchArea_Handle
try
nOut=nargout(MatchNURBSpatchArea_Handle);
catch
include_NURBSEngine
nOut=nargout(MatchNURBSpatchArea_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MatchNURBSpatchArea_Handle(varargin{:});
end
