function [varargout]=MakeCommentsFile(varargin)
% include_PostProcessing
global MakeCommentsFile_Handle
try
nOut=nargout(MakeCommentsFile_Handle);
catch
include_PostProcessing
nOut=nargout(MakeCommentsFile_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=MakeCommentsFile_Handle(varargin{:});
end
