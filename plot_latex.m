% Function turns the figure background white, increases font sizes, and
% allows mathematical formulae be inputted as figure titles. 
% (Makes figures appropriate for a Latex report)
function plot_latex(plt, xtitle, ytitle,ztitle, tl ,str)



%set the plot legend 
if isempty(str) == false 
    lgd = legend(str, 'Interpreter','tex','Orientation', 'Vertical'); 
    lgd.FontSize = 12;
end 


% set x label, y label and title of the graph, as interpreted using tex 
xlabel(xtitle,'fontsize',12,'fontweight','bold', 'Interpreter','tex');
ylabel(ytitle,'fontsize',12,'fontweight','bold','Interpreter','tex');
zlabel(ztitle,'fontsize',12,'fontweight','bold','Interpreter','tex');
title(tl,'fontsize',12,'fontweight','bold','Interpreter','tex');

% sets the fontsize of the values on the axes.
set(gca,'fontsize',12)

% the default matlab plot background is grey, this turns it to white 
set(gcf,'color','w');


% increase the width of lines plotted on the graph
% set(plt , 'linewidth',2);
set(findall(gca, 'Type', 'Line'),'LineWidth',2);

box on
grid minor

%fullscreen 
%figure('units','normalized','outerposition',[0 0 6 1])
fig=gcf;
fig.Position = [200, 200, 800,550];
end 
