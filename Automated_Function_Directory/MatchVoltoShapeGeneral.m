function [varargout]=MatchVoltoShapeGeneral(varargin)
% include_Optimisation
global MatchVoltoShapeGeneral_Handle
nOut=nargout(MatchVoltoShapeGeneral_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MatchVoltoShapeGeneral_Handle(varargin{:});
end
