% Written 10/24/03


% TO cleanup TDT system during debugging 

SBset([7 7],[0 0]);
RPhalt(RP);
RPclear(RP);
PAset(PAattns-PAattns+120);
SBset([],[]);
