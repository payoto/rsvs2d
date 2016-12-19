function [varargout]=InverseDesignBulk(varargin)
global InverseDesignBulk_Handle
nOut=nargout(InverseDesignBulk_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=InverseDesignBulk_Handle(varargin{:});
end
