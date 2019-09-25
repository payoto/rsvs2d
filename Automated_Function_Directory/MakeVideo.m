function [varargout]=MakeVideo(varargin)
% include_PostProcessing
global MakeVideo_Handle
try
nOut=nargout(MakeVideo_Handle);
catch
include_PostProcessing
nOut=nargout(MakeVideo_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeVideo_Handle(varargin{:});
end
