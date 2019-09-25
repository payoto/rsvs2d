function [varargout]=NacaOuterLimit4d(varargin)
% include_Optimisation
global NacaOuterLimit4d_Handle
try
nOut=nargout(NacaOuterLimit4d_Handle);
catch
include_Optimisation
nOut=nargout(NacaOuterLimit4d_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NacaOuterLimit4d_Handle(varargin{:});
end
