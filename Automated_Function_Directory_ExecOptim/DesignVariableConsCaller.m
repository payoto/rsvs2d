function [varargout]=DesignVariableConsCaller(varargin)
global DesignVariableConsCaller_Handle
nOut=nargout(DesignVariableConsCaller_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=DesignVariableConsCaller_Handle(varargin{:});
end
