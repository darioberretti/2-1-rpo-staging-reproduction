function [value, isterminal, direction] = event_xaxiscrossing(t, xx, isTerminal, evtDirection)
%EVENT_XAXISCROSSING Detect y = 0 crossings for the shared propagator.
%
% PROTOTYPE:
%   [value, isterminal, direction] = event_xaxiscrossing(t, xx, ...
%       isTerminal, evtDirection)
%
% DESCRIPTION:
%   Return the rotating-frame y coordinate as an ODE event value. For
%   t <= 1e-10 TU return 1 instead, suppressing the zero at the initial
%   epoch. The supplied time is absolute within the current integration.
%
% INPUT:
%   t             Scalar integration time [TU].
%   xx            State with y [DU] in its second component.
%   isTerminal    Logical controlling termination at a detected crossing.
%   evtDirection  Crossing direction: -1 (decreasing), 0 (either), 1 (increasing).
%
% OUTPUT:
%   value         Event value; zero identifies an accepted y = 0 crossing.
%   isterminal    Supplied terminal-event flag.
%   direction     Supplied crossing direction.

    if t > 1e-10
        value = xx(2);
    else
        value = 1;   
    end
    isterminal = isTerminal;
    direction = evtDirection;
end
