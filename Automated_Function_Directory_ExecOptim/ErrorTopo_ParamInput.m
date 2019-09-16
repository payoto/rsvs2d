function [varargout]=ErrorTopo_ParamInput(varargin)
% InverseDesign_ErrorTopo
global ErrorTopo_ParamInput_Handle
nOut=nargout(ErrorTopo_ParamInput_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=ErrorTopo_ParamInput_Handle(varargin{:});
end
