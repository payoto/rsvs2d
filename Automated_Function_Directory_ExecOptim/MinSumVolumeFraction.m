function [varargout]=MinSumVolumeFraction(varargin)
global MinSumVolumeFraction_Handle
nOut=nargout(MinSumVolumeFraction_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MinSumVolumeFraction_Handle(varargin{:});
end
