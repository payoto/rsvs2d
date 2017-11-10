function [varargout]=ExtractVolInfo(varargin)
global ExtractVolInfo_Handle
nOut=nargout(ExtractVolInfo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractVolInfo_Handle(varargin{:});
end
