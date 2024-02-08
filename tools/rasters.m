%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    % Routine to add rasters to current plot %
%    % data is an array of spike times        %
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rasters(data,color,event_draw,p)

if isempty(data)
  return
end

event_draw(:,p-2)=nan;
yl = get(gca,'ylim');
xl = get(gca,'xlim');

yt = get(gca,'ytick');
ytl = get(gca,'yticklabel');
NT = length(data(:,1));

unit  = yl(2)-yl(1);
offset = yl(2);
scale = 0.01*unit;
hold on
% line(xl,offset*[1 1],'color','k')
for n=1:NT
  indx = find(data(n,:)>xl(1) & data(n,:)<xl(2) ...
                              & data(n,:) ~= 0);
  yloc = offset + n*scale;
  if ~isempty(indx)
    plot(data(n,indx),yloc*ones(length(indx),1),[color,'.'],'MarkerSize',1.5);
    hold on;
    plot(event_draw(n,1),yloc,['m','.'],'MarkerSize',4); hold on;
    plot(event_draw(n,2),yloc,['g','.'],'MarkerSize',4); hold on;
    if size(event_draw,2)>3
    plot(event_draw(n,4),yloc,['c','.'],'MarkerSize',4); hold on;
    plot(event_draw(n,5),yloc,['y','.'],'MarkerSize',4); hold on;
    end  
  end
end

set(gca,'ylim',[yl(1) yloc+scale]);  
set(gca,'ytick',yt);
set(gca,'yticklabel',ytl);
plot([0 0],get(gca,'ylim'),'k');
hold off 








