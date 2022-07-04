function []  = Plot_cows(health_cows_local,infected_cows_local,size_farm,t)

% this is the function to plot the location of cows


figure
x1=health_cows_local(:,1);
y1=health_cows_local(:,2);
x2=infected_cows_local(:,1);
y2=infected_cows_local(:,2);
leng=size_farm(1,2);
axis square
hold on
plot(x2,y2,'x','color','r')   
plot(x1,y1,'.','color','b')     
text(-leng*0.775,-leng*0.625,'drinking area','FontSize',8)
text(-leng*0.775,leng*0.625,'drinking area','FontSize',8)
text(-leng*0.125,leng*0.05,'drinking area','FontSize',8)
text(leng*0.475,leng*0.625,'drinking area','FontSize',8)
text(leng*0.475,leng*-0.625,'drinking area','FontSize',8)
rectangle('position',[-leng*0.75 -leng*0.75 leng/4 leng/4],'LineWidth',2,'LineStyle',':');  
rectangle('position',[-leng*0.75 leng*0.5 leng/4 leng/4],'LineWidth',2,'LineStyle',':');  
rectangle('position',[leng*0.5 leng*0.5 leng/4 leng/4],'LineWidth',2,'LineStyle',':');  
rectangle('position',[-leng*0.1 -leng*0.1 leng/4 leng/4],'LineWidth',2,'LineStyle',':');  
rectangle('position',[leng*0.5 -leng*0.75 leng/4 leng/4],'LineWidth',2,'LineStyle',':');   
set(gca,'XLim',[size_farm(1,1) size_farm(1,2)]);        
set(gca,'XTick',[size_farm(1,1):size_farm(1,2)/10:size_farm(1,2)]);  
set(gca,'YLim',[size_farm(1,1) size_farm(1,2)]);        
set(gca,'YTick',[size_farm(1,1):size_farm(1,2)/10:size_farm(1,2)]);
if t==0
   title(' Initial location of healthy and infected cows ');
else
   title(['the location of cows in day ',num2str(t/24)]);
end
xlabel('represent 100m for each num ')
ylabel('represent 100m for each num ')


legend('infected cows','Health cows')

end

