function [varargout]=ConstantArea_Klunker(varargin)
global ConstantArea_Klunker_Handle
nOut=nargout(ConstantArea_Klunker_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConstantArea_Klunker_Handle(varargin{:});
end
