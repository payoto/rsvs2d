function [varargout]=ExtractDesignVariableConnectivity(varargin)
global ExtractDesignVariableConnectivity_Handle
nOut=nargout(ExtractDesignVariableConnectivity_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ExtractDesignVariableConnectivity_Handle(varargin{:});
end
