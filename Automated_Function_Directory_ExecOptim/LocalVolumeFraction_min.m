function [varargout]=LocalVolumeFraction_min(varargin)
global LocalVolumeFraction_min_Handle
nOut=nargout(LocalVolumeFraction_min_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LocalVolumeFraction_min_Handle(varargin{:});
end
