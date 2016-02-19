function [varargout]=ModifReshape(varargin)
global ModifReshape_Handle
nOut=nargout(ModifReshape_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ModifReshape_Handle(varargin{:});
end
