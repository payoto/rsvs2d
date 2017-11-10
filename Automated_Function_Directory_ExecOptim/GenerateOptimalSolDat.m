function [varargout]=GenerateOptimalSolDat(varargin)
global GenerateOptimalSolDat_Handle
nOut=nargout(GenerateOptimalSolDat_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateOptimalSolDat_Handle(varargin{:});
end
