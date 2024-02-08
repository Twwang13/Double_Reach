function Fig=Exam_Cell
set(0,'DefaultAxesFontName','Arial');
set(0,'defaultTextFontName','Arial');
set(0,'defaultUipanelFontName','Arial');

 


[filename, pathname] = uigetfile( {'*.*'}, 'Pick a file','on');
[spike_processed,event] = process_data(pathname,filename); % replace with function D_process_data for processed43_unit2.mat
event(:,4) = event(:,4)-0.15; % correction for touch screen delay 150ms
zone = 0.5;
plot_name = [];
tt_analysis_plot(spike_processed,event,plot_name,zone);

%% 
h = get(gcf);
ax6 = h.Children(6);
anchor1 = get( ax6,'Position');
% set(ax6,'position',[anchor1(1),anchor1(2),anchor1(3),anchor1(4)]);

x = 0.1221;
hght = 0.3;

ax7 = h.Children(7);
po = get( ax7, 'Position' );
set(ax7,'position',[anchor1(1)+x,anchor1(2)+hght,po(3),po(4)]);

ax8 = h.Children(8);
po = get( ax8, 'Position' );
set(ax8,'position',[anchor1(1)+2.7*x,anchor1(2)+hght,po(3),po(4)]);

ax4 = h.Children(4);
po = get( ax4, 'Position' );
set(ax4,'position',[anchor1(1)+2.7*x,anchor1(2)-hght,po(3),po(4)]);

ax5 = h.Children(5);
po = get( ax5, 'Position' );
set(ax5,'position',[anchor1(1)+x,anchor1(2)-hght,po(3),po(4)]);

ax9 = h.Children(9);
po = get(ax9,'Position');
set(ax9,'position',[anchor1(1)+3.7*x,anchor1(2),po(3),po(4)]);

Fig = ['3GO ','4MO ','5ME '];
end