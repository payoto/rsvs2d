function [varargout]=TestNormalOutPointing(varargin)
% include_NURBSEngine
global TestNormalOutPointing_Handle
try
nOut=nargout(TestNormalOutPointing_Handle);
catch
include_NURBSEngine
nOut=nargout(TestNormalOutPointing_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=TestNormalOutPointing_Handle(varargin{:});
end
