function [varargout]=OpenTecLayFile(varargin)
% include_PostProcessing
global OpenTecLayFile_Handle
nOut=nargout(OpenTecLayFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenTecLayFile_Handle(varargin{:});
end
