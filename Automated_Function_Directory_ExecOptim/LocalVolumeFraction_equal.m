function [varargout]=LocalVolumeFraction_equal(varargin)
global LocalVolumeFraction_equal_Handle
nOut=nargout(LocalVolumeFraction_equal_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=LocalVolumeFraction_equal_Handle(varargin{:});
end
