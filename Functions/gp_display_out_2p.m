function astroGpe = gp_display_out_2p(Gpe,labeled,scaleBar)
% Displays the groupe of astrocytes specified in 'Gpe'
% 25/03/2019
astroGpe = zeros(size(labeled));
for Ag = 1:length(Gpe)
    Agrp = Gpe(Ag);
%         x = find(labeled == Ag); 
    astroGpe(labeled == Agrp)=1;
end    
% Add scale bar 100µm
% insertShape(astroGpe,'Line',[20 20,100 20],'LineWidth',5);
if strcmp(scaleBar,'on')
    scaleBar = zeros(size(astroGpe));
    scaleBar(20:22,20:100) = 2;
    astroGpe = astroGpe + scaleBar;
end
end
