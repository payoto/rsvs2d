function [varargout]=ImposePresets(varargin)
global ImposePresets_Handle
nOut=nargout(ImposePresets_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ImposePresets_Handle(varargin{:});
end
