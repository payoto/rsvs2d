function [varargout]=RobustSymmetryAssignement(varargin)
global RobustSymmetryAssignement_Handle
nOut=nargout(RobustSymmetryAssignement_Handle);
nOutReq=nargout;
nOut(nOut<0)=nOutReq;
[varargout{1:nOut}]=RobustSymmetryAssignement_Handle(varargin{:});
end
