function [varargout]=ConvertProfToFill(varargin)
global ConvertProfToFill_Handle
nOut=nargout(ConvertProfToFill_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ConvertProfToFill_Handle(varargin{:});
end
