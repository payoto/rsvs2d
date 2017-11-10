function [varargout]=GenerateIterResultBinary(varargin)
% OptimisationOutput
global GenerateIterResultBinary_Handle
nOut=nargout(GenerateIterResultBinary_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateIterResultBinary_Handle(varargin{:});
end
