function [varargout]=GenerateValidationData(varargin)
% include_Validation
global GenerateValidationData_Handle
try
nOut=nargout(GenerateValidationData_Handle);
catch
include_Validation
nOut=nargout(GenerateValidationData_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateValidationData_Handle(varargin{:});
end
