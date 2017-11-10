function [varargout]=OpenParamFile(varargin)
% include_PostProcessing
global OpenParamFile_Handle
nOut=nargout(OpenParamFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OpenParamFile_Handle(varargin{:});
end
