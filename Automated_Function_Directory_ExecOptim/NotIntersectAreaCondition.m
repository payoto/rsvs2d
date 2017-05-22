function [varargout]=NotIntersectAreaCondition(varargin)
global NotIntersectAreaCondition_Handle
nOut=nargout(NotIntersectAreaCondition_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NotIntersectAreaCondition_Handle(varargin{:});
end
