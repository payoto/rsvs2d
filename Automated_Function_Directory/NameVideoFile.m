function [varargout]=NameVideoFile(varargin)
% include_PostProcessing
global NameVideoFile_Handle
nOut=nargout(NameVideoFile_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=NameVideoFile_Handle(varargin{:});
end
