function [varargout]=GenerateVariableString(varargin)
% include_PostProcessing
global GenerateVariableString_Handle
nOut=nargout(GenerateVariableString_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateVariableString_Handle(varargin{:});
end
