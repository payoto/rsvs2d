function [varargout]=MeanVolumeFraction(varargin)
global MeanVolumeFraction_Handle
nOut=nargout(MeanVolumeFraction_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MeanVolumeFraction_Handle(varargin{:});
end
