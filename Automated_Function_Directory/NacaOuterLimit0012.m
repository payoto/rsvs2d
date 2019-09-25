function [varargout]=NacaOuterLimit0012(varargin)
% include_Optimisation
global NacaOuterLimit0012_Handle
try
nOut=nargout(NacaOuterLimit0012_Handle);
catch
include_Optimisation
nOut=nargout(NacaOuterLimit0012_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NacaOuterLimit0012_Handle(varargin{:});
end
