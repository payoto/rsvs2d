function [varargout]=ExtractSnaxelLoops(varargin)
global ExtractSnaxelLoops_Handle
nOut=nargout(ExtractSnaxelLoops_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractSnaxelLoops_Handle(varargin{:});
end
