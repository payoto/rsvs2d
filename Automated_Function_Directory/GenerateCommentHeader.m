function [varargout]=GenerateCommentHeader(varargin)
% include_PostProcessing
global GenerateCommentHeader_Handle
nOut=nargout(GenerateCommentHeader_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateCommentHeader_Handle(varargin{:});
end
