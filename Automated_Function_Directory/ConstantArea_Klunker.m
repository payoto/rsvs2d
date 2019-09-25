function [varargout]=ConstantArea_Klunker(varargin)
% include_Optimisation
global ConstantArea_Klunker_Handle
try
nOut=nargout(ConstantArea_Klunker_Handle);
catch
include_Optimisation
nOut=nargout(ConstantArea_Klunker_Handle);
end
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConstantArea_Klunker_Handle(varargin{:});
end
