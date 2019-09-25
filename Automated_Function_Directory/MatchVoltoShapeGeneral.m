function [varargout]=MatchVoltoShapeGeneral(varargin)
% include_Optimisation
global MatchVoltoShapeGeneral_Handle
try
nOut=nargout(MatchVoltoShapeGeneral_Handle);
catch
include_Optimisation
nOut=nargout(MatchVoltoShapeGeneral_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MatchVoltoShapeGeneral_Handle(varargin{:});
end
