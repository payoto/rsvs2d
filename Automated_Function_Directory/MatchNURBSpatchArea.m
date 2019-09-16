function [varargout]=MatchNURBSpatchArea(varargin)
% include_NURBSEngine
global MatchNURBSpatchArea_Handle
nOut=nargout(MatchNURBSpatchArea_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MatchNURBSpatchArea_Handle(varargin{:});
end
