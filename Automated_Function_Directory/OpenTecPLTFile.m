function [varargout]=OpenTecPLTFile(varargin)
global OpenTecPLTFile_Handle
nOut=nargout(OpenTecPLTFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenTecPLTFile_Handle(varargin{:});
end
