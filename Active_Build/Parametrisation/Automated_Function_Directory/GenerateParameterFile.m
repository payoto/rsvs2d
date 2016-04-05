function [varargout]=GenerateParameterFile(varargin)
global GenerateParameterFile_Handle
nOut=nargout(GenerateParameterFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateParameterFile_Handle(varargin{:});
end
