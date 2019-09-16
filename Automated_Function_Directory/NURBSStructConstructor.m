function [varargout]=NURBSStructConstructor(varargin)
% include_NURBSEngine
global NURBSStructConstructor_Handle
nOut=nargout(NURBSStructConstructor_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NURBSStructConstructor_Handle(varargin{:});
end
