function [varargout]=ProjectColormap(varargin)
global ProjectColormap_Handle
nOut=nargout(ProjectColormap_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ProjectColormap_Handle(varargin{:});
end
