function [varargout]=GenerateProfileBinary(varargin)
% include_PostProcessing
global GenerateProfileBinary_Handle
try
nOut=nargout(GenerateProfileBinary_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateProfileBinary_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateProfileBinary_Handle(varargin{:});
end
