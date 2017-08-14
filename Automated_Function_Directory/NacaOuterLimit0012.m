function [varargout]=NacaOuterLimit0012(varargin)
% include_Optimisation
global NacaOuterLimit0012_Handle
nOut=nargout(NacaOuterLimit0012_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NacaOuterLimit0012_Handle(varargin{:});
end
