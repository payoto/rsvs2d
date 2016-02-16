function [varargout]=ModifReshape(varargin)
global ModifReshape_Handle
nOut=nargout(ModifReshape_Handle);
[varargout{1:nOut}]=ModifReshape_Handle(varargin{:});
end
