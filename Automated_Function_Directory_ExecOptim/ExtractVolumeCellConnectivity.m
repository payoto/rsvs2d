function [varargout]=ExtractVolumeCellConnectivity(varargin)
global ExtractVolumeCellConnectivity_Handle
nOut=nargout(ExtractVolumeCellConnectivity_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVolumeCellConnectivity_Handle(varargin{:});
end
