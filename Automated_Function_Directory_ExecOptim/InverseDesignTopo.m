function [varargout]=InverseDesignTopo(varargin)
global InverseDesignTopo_Handle
nOut=nargout(InverseDesignTopo_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InverseDesignTopo_Handle(varargin{:});
end
