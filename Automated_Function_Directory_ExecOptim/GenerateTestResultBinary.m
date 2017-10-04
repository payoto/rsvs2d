function [varargout]=GenerateTestResultBinary(varargin)
% OptimisationOutput
global GenerateTestResultBinary_Handle
nOut=nargout(GenerateTestResultBinary_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateTestResultBinary_Handle(varargin{:});
end
