function [varargout]=WriteFullInfoProfile(varargin)
% include_PostProcessing
global WriteFullInfoProfile_Handle
try
nOut=nargout(WriteFullInfoProfile_Handle);
catch
include_PostProcessing
nOut=nargout(WriteFullInfoProfile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=WriteFullInfoProfile_Handle(varargin{:});
end
