function [varargout]=GenerateVariableString(varargin)
% include_PostProcessing
global GenerateVariableString_Handle
try
nOut=nargout(GenerateVariableString_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateVariableString_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateVariableString_Handle(varargin{:});
end
