function [SCCpos,SCCs_belowBF,SCCs_aboveBF,BFind,centerBF_kHz]=getSCCindsWithBF(unitBF_kHz,NSAC_BFs_kHz,NSCC_BFs_kHz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracts from NSCC all SCCinds that have unitBF_kHz as one of the freq
% pairs. If unitBF_kHz is not part of the NSAC inventory (NSAC_BFs_kHz), centerBF_kHz 
% substitutes unitBF_kHz with the closest freq found in the inventory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NSAC_BFs_kHz_vec=cell2mat(NSAC_BFs_kHz');
[y,BFind]=min(abs(NSAC_BFs_kHz_vec-unitBF_kHz));

centerBF_kHz=NSAC_BFs_kHz_vec(BFind);
[SCCrow,SCCcol]=find(cell2mat(NSCC_BFs_kHz')==centerBF_kHz);
SCCpos=sortrows([SCCrow,SCCcol]);%SCCcol entries are either 1 or 2.
SCCs_belowBF=(SCCpos(find(SCCpos(:,2)==2),1));
SCCs_aboveBF=(SCCpos(find(SCCpos(:,2)==1),1));
return;