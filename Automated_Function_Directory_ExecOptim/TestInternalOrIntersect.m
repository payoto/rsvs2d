function [varargout]=TestInternalOrIntersect(varargin)
global TestInternalOrIntersect_Handle
nOut=nargout(TestInternalOrIntersect_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TestInternalOrIntersect_Handle(varargin{:});
end
