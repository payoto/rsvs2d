function [varargout]=ImageProcess(varargin)
% include_Utilities
global ImageProcess_Handle
nOut=nargout(ImageProcess_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ImageProcess_Handle(varargin{:});
end
