function [varargout]=FindProfile(varargin)
% include_PostProcessing
global FindProfile_Handle
try
nOut=nargout(FindProfile_Handle);
catch
include_PostProcessing
nOut=nargout(FindProfile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FindProfile_Handle(varargin{:});
end
