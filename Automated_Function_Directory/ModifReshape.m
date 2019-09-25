function [varargout]=ModifReshape(varargin)
% include_SnakeParam
global ModifReshape_Handle
try
nOut=nargout(ModifReshape_Handle);
catch
include_SnakeParam
nOut=nargout(ModifReshape_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ModifReshape_Handle(varargin{:});
end
