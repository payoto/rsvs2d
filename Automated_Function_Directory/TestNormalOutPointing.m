function [varargout]=TestNormalOutPointing(varargin)
% include_NURBSEngine
global TestNormalOutPointing_Handle
nOut=nargout(TestNormalOutPointing_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TestNormalOutPointing_Handle(varargin{:});
end
