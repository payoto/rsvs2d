function [varargout]=GenerateValidationData(varargin)
global GenerateValidationData_Handle
nOut=nargout(GenerateValidationData_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateValidationData_Handle(varargin{:});
end
