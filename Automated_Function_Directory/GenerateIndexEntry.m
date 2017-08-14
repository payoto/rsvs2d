function [varargout]=GenerateIndexEntry(varargin)
% include_PostProcessing
global GenerateIndexEntry_Handle
nOut=nargout(GenerateIndexEntry_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateIndexEntry_Handle(varargin{:});
end
