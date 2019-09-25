function [varargout]=GenerateParameterFile(varargin)
% include_PostProcessing
global GenerateParameterFile_Handle
try
nOut=nargout(GenerateParameterFile_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateParameterFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateParameterFile_Handle(varargin{:});
end
