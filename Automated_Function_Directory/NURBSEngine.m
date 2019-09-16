function [varargout]=NURBSEngine(varargin)
% include_NURBSEngine
global NURBSEngine_Handle
nOut=nargout(NURBSEngine_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NURBSEngine_Handle(varargin{:});
end
