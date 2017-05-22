function [varargout]=NotIntersectCondition(varargin)
global NotIntersectCondition_Handle
nOut=nargout(NotIntersectCondition_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NotIntersectCondition_Handle(varargin{:});
end
