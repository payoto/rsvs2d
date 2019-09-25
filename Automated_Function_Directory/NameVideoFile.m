function [varargout]=NameVideoFile(varargin)
% include_PostProcessing
global NameVideoFile_Handle
try
nOut=nargout(NameVideoFile_Handle);
catch
include_PostProcessing
nOut=nargout(NameVideoFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NameVideoFile_Handle(varargin{:});
end
