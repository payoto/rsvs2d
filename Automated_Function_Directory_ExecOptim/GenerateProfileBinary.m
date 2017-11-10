function [varargout]=GenerateProfileBinary(varargin)
% OptimisationOutput
global GenerateProfileBinary_Handle
nOut=nargout(GenerateProfileBinary_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateProfileBinary_Handle(varargin{:});
end
