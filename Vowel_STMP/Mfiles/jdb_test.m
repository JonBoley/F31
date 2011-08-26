%%
SMP_NSCC_Oct_at_Rho=cell(length(FeatINDs),length(Nattens_dB));
SMP_NSCC_CD_at_Oct=cell(length(FeatINDs),length(Nattens_dB));


for ATTind=1:length(Nattens_dB)
    ROWind=0;
    for FeatIND=FeatINDs
        ROWind=ROWind+1;
        SMP_NSCC_Oct_at_Rho{ROWind,ATTind}=NaN;
        SMP_NSCC_CD_at_Oct{ROWind,ATTind}=NaN;
        
        [SCCpos,SCCs_belowBF,SCCs_aboveBF,BFind,BFs(ROWind)]=getSCCindsWithBF(unit.Info.BF_kHz,NSAC_BFs_kHz{ROWind,ATTind},NSCC_BFs_kHz{ROWind,ATTind});
        for SCCposind=1:size(SCCpos,1)
            if ~isnan(NSCC_Rho{ROWind,ATTind}{SCCpos(SCCposind,1)})
                MyNSCC(SCCposind) = NSCC_Rho{ROWind,ATTind}{SCCpos(SCCposind,1)};
                Octaves(SCCposind) = log2(NSCC_BFs_kHz{ROWind,ATTind}{SCCpos(SCCposind,1)}(3-SCCpos(SCCposind,2))/BFs(ROWind));
                MyCD(SCCposind) = NSCC_CDs_usec{ROWind,ATTind}{SCCpos(SCCposind,1)};
            end
        end
        Octaves_interp = min(Octaves):0.01:max(Octaves);
        MyNSCC_interp = interp1(Octaves,MyNSCC,Octaves_interp);
        MyCD_interp = interp1(Octaves,MyCD,Octaves_interp);
        
        % find CD at 0.25 octaves away from feature
        SMP_NSCC_CD_at_Oct{ROWind,ATTind} = max(abs(MyCD_interp(abs(Octaves_interp)<=0.25)));
        % find oct diff at which NSCC_Rho=0.6
        SMP_NSCC_Oct_at_Rho{ROWind,ATTind} = max(abs(Octaves_interp(MyNSCC_interp>=0.6)));
        if isempty(SMP_NSCC_Oct_at_Rho{ROWind,ATTind})
            SMP_NSCC_Oct_at_Rho{ROWind,ATTind} = max(abs(Octaves));
        end
        clear MyNSCC MyNSCC_interp Octaves Octaves_interp MyCD;
    end
end
