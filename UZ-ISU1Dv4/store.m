% Part of a new version of the solution for the Richards equation
% Written May 2012 by Daniel Erdal
% Compacted version written May 2012

% Version 1.0
% Last Documented edition: 30.05.12 (DE)

function storeX=store(h,wc_new,wc_old,param,dt,delzCS,head_old)


% Calculate the storage and the first term for the evaluation of the
% richards equation.

storeX=(wc_new./param.poro.*param.storage.*(h-head_old)+wc_new-wc_old).*delzCS./dt;

