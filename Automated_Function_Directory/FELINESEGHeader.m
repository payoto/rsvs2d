function [varargout]=FELINESEGHeader(varargin)
% include_PostProcessing
global FELINESEGHeader_Handle
try
nOut=nargout(FELINESEGHeader_Handle);
catch
include_PostProcessing
nOut=nargout(FELINESEGHeader_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=FELINESEGHeader_Handle(varargin{:});
end
