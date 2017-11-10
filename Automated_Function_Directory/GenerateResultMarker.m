function [varargout]=GenerateResultMarker(varargin)
% include_PostProcessing
global GenerateResultMarker_Handle
nOut=nargout(GenerateResultMarker_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=GenerateResultMarker_Handle(varargin{:});
end
