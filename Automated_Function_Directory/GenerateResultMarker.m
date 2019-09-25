function [varargout]=GenerateResultMarker(varargin)
% include_PostProcessing
global GenerateResultMarker_Handle
try
nOut=nargout(GenerateResultMarker_Handle);
catch
include_PostProcessing
nOut=nargout(GenerateResultMarker_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateResultMarker_Handle(varargin{:});
end
