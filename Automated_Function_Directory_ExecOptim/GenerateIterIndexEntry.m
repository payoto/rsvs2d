function [varargout]=GenerateIterIndexEntry(varargin)
% OptimisationOutput
global GenerateIterIndexEntry_Handle
nOut=nargout(GenerateIterIndexEntry_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateIterIndexEntry_Handle(varargin{:});
end
