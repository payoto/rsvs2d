function [varargout]=GenerateRestartBinary(varargin)
% include_PostProcessing
global GenerateRestartBinary_Handle
nOut=nargout(GenerateRestartBinary_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateRestartBinary_Handle(varargin{:});
end
