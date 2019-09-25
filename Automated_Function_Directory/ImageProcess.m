function [varargout]=ImageProcess(varargin)
% include_Utilities
global ImageProcess_Handle
try
nOut=nargout(ImageProcess_Handle);
catch
include_Utilities
nOut=nargout(ImageProcess_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ImageProcess_Handle(varargin{:});
end
