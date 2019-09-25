function [varargout]=Image2Fill(varargin)
% include_Utilities
global Image2Fill_Handle
try
nOut=nargout(Image2Fill_Handle);
catch
include_Utilities
nOut=nargout(Image2Fill_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=Image2Fill_Handle(varargin{:});
end
