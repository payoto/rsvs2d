function [varargout]=ExtractSnakInfo_vertflow(varargin)
% include_Optimisation
global ExtractSnakInfo_vertflow_Handle
nOut=nargout(ExtractSnakInfo_vertflow_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractSnakInfo_vertflow_Handle(varargin{:});
end
