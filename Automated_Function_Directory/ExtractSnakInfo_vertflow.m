function [varargout]=ExtractSnakInfo_vertflow(varargin)
% include_Optimisation
global ExtractSnakInfo_vertflow_Handle
try
nOut=nargout(ExtractSnakInfo_vertflow_Handle);
catch
include_Optimisation
nOut=nargout(ExtractSnakInfo_vertflow_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractSnakInfo_vertflow_Handle(varargin{:});
end
