function [varargout]=GenerateIndexEntry(varargin)
% include_PostProcessing
global GenerateIndexEntry_Handle
try
nOut=nargout(GenerateIndexEntry_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateIndexEntry_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateIndexEntry_Handle(varargin{:});
end
