function [varargout]=GenerateRestartBinary(varargin)
% include_PostProcessing
global GenerateRestartBinary_Handle
try
nOut=nargout(GenerateRestartBinary_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateRestartBinary_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateRestartBinary_Handle(varargin{:});
end
