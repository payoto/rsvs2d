function [varargout]=GenerateOptimalSolDir(varargin)
% OptimisationOutput
global GenerateOptimalSolDir_Handle
nOut=nargout(GenerateOptimalSolDir_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateOptimalSolDir_Handle(varargin{:});
end
