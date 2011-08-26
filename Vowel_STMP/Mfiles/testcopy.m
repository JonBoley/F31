%%
hfigure_orig=figure(1000);
haxes_orig = findobj(hfigure_orig, 'Type', 'Axes');
haxes_children = allchild(haxes_orig);%, 'Children');

figure(501), set(gcf,'Name','Comparison of Neural Coding Strategies'); 
clf;
set(gcf,'units','norm','pos',[0.2234    0.7119-0.2344    0.4297    0.2344])%,'Resize','off')
% copyobj(haxes_orig, gcf); % copy everything over to this new figure

% Spatio-temporal plots
copyobj(haxes_children{end}, subplot(2,5,1));
xlim(XLIMITS_perhist)
set(gca,'YTick',YTICKS,'YTickLabel',YTICKS)
ylim(YLIMITS)  % Same Ylimits for all plots
                
% Rate plots
copyobj(haxes_children{end-2}, subplot(2,5,2));
xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]\nO: ALSR'),'FontSize',6)
xlim(XLIMITS_rate)
set(gca,'XDir','reverse')
set(gca,'YTick',YTICKS,'YTickLabel',YTICKS)
ylim(YLIMITS)  % Same Ylimits for all plots

copyobj(haxes_children{end-3}, subplot(2,5,3));
% subplot(2,5,4) gets rho_tfs for F1
% subplot(2,5,5) gets CD for F1
copyobj(haxes_children{end-5}, subplot(2,5,6));
copyobj(haxes_children{end-6}, subplot(2,5,7));
copyobj(haxes_children{end-7}, subplot(2,5,8));
% subplot(2,5,9) gets rho_tfs for T1
% subplot(2,5,10) gets CD for T1

haxes = findobj(gcf, 'Type', 'Axes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% copied from other file
NUMrows=2;
Xcorner=0.05;
Xwidth1=.5;
Xshift1=0.05;
Xwidth2=.1;
Xshift2=0.03;
Ycorner=0.05;
Yshift=0.07;
Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2
TICKlength=0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(haxes(end),'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
set(haxes(end-1),'pos',[0.4451    0.8913    0.0942    0.0473]) % legend
set(haxes(end-2),'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(haxes(end-3),'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(haxes(end-4),'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(haxes(end-5),'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
set(haxes(end-6),'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(haxes(end-7),'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(haxes(end-8),'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
