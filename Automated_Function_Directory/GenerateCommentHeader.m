function [varargout]=GenerateCommentHeader(varargin)
% include_PostProcessing
global GenerateCommentHeader_Handle
try
nOut=nargout(GenerateCommentHeader_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateCommentHeader_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateCommentHeader_Handle(varargin{:});
end
