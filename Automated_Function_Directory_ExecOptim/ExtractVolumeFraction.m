function [varargout]=ExtractVolumeFraction(varargin)
global ExtractVolumeFraction_Handle
nOut=nargout(ExtractVolumeFraction_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVolumeFraction_Handle(varargin{:});
end
