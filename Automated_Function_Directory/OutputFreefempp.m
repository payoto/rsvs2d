function [varargout]=OutputFreefempp(varargin)
% include_PostProcessing
global OutputFreefempp_Handle
try
nOut=nargout(OutputFreefempp_Handle);
catch
include_PostProcessing
nOut=nargout(OutputFreefempp_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=OutputFreefempp_Handle(varargin{:});
end
