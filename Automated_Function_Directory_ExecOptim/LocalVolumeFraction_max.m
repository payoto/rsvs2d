function [varargout]=LocalVolumeFraction_max(varargin)
global LocalVolumeFraction_max_Handle
nOut=nargout(LocalVolumeFraction_max_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LocalVolumeFraction_max_Handle(varargin{:});
end
