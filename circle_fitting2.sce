clear
clearglobal
//clc
global polygon
global circle
function draw_circle(c)
    theta=linspace(0,2*%pi,16);
    x=c.pos(1)+c.radius*cos(theta);
    y=c.pos(2)+c.radius*sin(theta);
    xstring(c.pos(1),c.pos(2),string(i));
    fa=0.15;
    fv=0.1;
    xarrows([c.pos(1) c.pos(1)+fa*c.accel(1)] ,[c.pos(2) c.pos(2)+fa*c.accel(2)])
    xarrows([c.pos(1) c.pos(1)+fv*c.vel(1)] ,[c.pos(2) c.pos(2)+fv*c.vel(2)],-1,3)
    plot(x,y);
endfunction
function findAllDistances()
    global polygon
    global circle
    for i=1:size(circle,1)
        circle(i)=findCircDistances(circle(i));
//        linedist=[];
//        circdist=1e5*ones(size(circle,1),size(circle,1));
//        x0=circle(i).pos(1);y0=circle(i).pos(2);
//        totangle=0;
//        for j=1:size(polygon.lines,1)
//            x1=polygon.lines(j,1);y1=polygon.lines(j,2);
//            x2=polygon.lines(j,3);y2=polygon.lines(j,4);
//            linedist=[linedist abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2)];
//            anglediff=atand(y1-y0,x1-x0) - atand(y2-y0,x2-x0);
//            totangle=totangle + anglediff;
//            if anglediff > 180
//                totangle=totangle-360;
//            elseif anglediff < -180
//                totangle=totangle+360;
//            end
//        end
//        circle(i).linedist=linedist;
//        if totangle<0.01 | min(linedist) <= circle(i).radius
//            circle(i).inside_polygon=0;
//        else
//            circle(i).inside_polygon=1;
//        end
//        for j=i+1:size(circle,1)
//            circdist(i,j)=[sqrt((circle(i).pos(1)-circle(j).pos(1))^2 + (circle(i).pos(2)-circle(j).pos(2))^2)];
//            circdist(j,i)=circdist(i,j);
//        end
//        circle(i).circdist=circdist;
//        if  min(circdist) <= 2*circle(i).radius
//            circle(i).clearcircles=0;
//        else
//            circle(i).clearcircles=1;
//        end
    end
endfunction
function c=findCircDistances(cold)
    global polygon
    global circle
    c=cold;
    linedist=[];
    circdist=zeros(size(circle,1),1);
    x0=c.pos(1);y0=c.pos(2);
    totangle=0;
    for j=1:size(polygon.lines,1)
        x1=polygon.lines(j,1);y1=polygon.lines(j,2);
        x2=polygon.lines(j,3);y2=polygon.lines(j,4);
        linedist=[linedist abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2)];
        anglediff= - atand(y1-y0,x1-x0) + atand(y2-y0,x2-x0);
        totangle=totangle + anglediff;
        if anglediff > 180
            totangle=totangle-360;
        elseif anglediff < -180
            totangle=totangle+360;
        end
    end
    c.linedist=linedist;
    //pause;
    if  min(linedist) < cold.radius// | totangle<0.01 
        c.inside_polygon=0;
    else
        c.inside_polygon=1;
    end
    for j=1:size(circle,1)
        circdist(j)=[sqrt((c.pos(1)-circle(j).pos(1))^2 + (c.pos(2)-circle(j).pos(2))^2)];
    end
    c.circdist=circdist;
    if  min(circdist(find(circdist>0))) < 2*c.radius
        c.clearofcircles=0; 
    else
        c.clearofcircles=1;
    end
endfunction
function c=dampingforce(c)
    k=0.6;pow=1;
    //c.accel=c.accel-k*(c.vel^pow)';

endfunction
function c=gravityforce(c)
    gravity.accel=[0 -9.8];
    c.accel=gravity.accel;
endfunction
function c=new_pos_vel(c,dt)
    c1=c;
    c1.vel=c1.vel+c1.accel'*dt;
    c1.pos=c1.pos+c1.vel'*dt;
    c=c1;
endfunction
function c=getnewspeed(cold,coldold,lineangle,elasticity)
    c=cold;
    try
        //pause
        ax=-9.8*sind(lineangle);
        ay=-9.8*cosd(lineangle);
        px=cold.pos(1)*cosd(lineangle)+cold.pos(2)*sind(lineangle);
        py=-cold.pos(1)*sind(lineangle)+cold.pos(2)*cosd(lineangle);
        //pxold=coldold.pos(1)*cosd(lineangle)+coldold.pos(2)*sind(lineangle);
        //pyold=-coldold.pos(1)*sind(lineangle)+coldold.pos(2)*cosd(lineangle);
        vx=cold.vel(1)*cosd(lineangle)+cold.vel(2)*sind(lineangle);
        vy=-cold.vel(1)*sind(lineangle)+cold.vel(2)*cosd(lineangle);
        vy=-elasticity*vy;
        
        //ay=-ay;

        c.vel(1)=vx*cosd(lineangle)-vy*sind(lineangle);
        c.vel(2)=vx*sind(lineangle)+vy*cosd(lineangle);
        c.pos(1)=px*cosd(lineangle)-py*sind(lineangle);
        c.pos(2)=px*sind(lineangle)+py*cosd(lineangle);
        c.pos(2)=coldold.pos(2);
        
        c.accel(1)=c.accel(1) - (ax*cosd(lineangle)-ay*sind(lineangle));
        c.accel(2)=c.accel(2) - (ax*sind(lineangle)+ay*cosd(lineangle));
//pause
    catch
        disp('in the line-circle function')
        pause
    end

endfunction
function [c1new,c2new]=getnewspeedcircles(c1,c1old,c2,elasticity)
    c1new=c1;
    c2new=c2;
    try
        normalangle=atand((c1.pos(2)-c2.pos(2))/(c1.pos(1)-c2.pos(1)));
        lineangle=90-normalangle;
        //pause;
        ax1=-9.8*sind(lineangle);
        ay1=-9.8*cosd(lineangle);
        vx1=c1.vel(1)*cosd(lineangle)-c1.vel(2)*sind(lineangle);
        vy1=c1.vel(1)*sind(lineangle)+c1.vel(2)*cosd(lineangle);
        vx2=c2.vel(1)*cosd(lineangle)-c2.vel(2)*sind(lineangle);
        vy2=c2.vel(1)*sind(lineangle)+c2.vel(2)*cosd(lineangle);
        
        tmp1=elasticity*vy1;
        tmp2=elasticity*vy2;
        vy1=tmp2;
        vy2=tmp1;

        c1new.vel(1)=vx1*cosd(lineangle)+vy1*sind(lineangle);
        c1new.vel(2)=-vx1*sind(lineangle)+vy1*cosd(lineangle);
        c2new.vel(1)=vx2*cosd(lineangle)+vy2*sind(lineangle);
        c2new.vel(2)=-vx2*sind(lineangle)+vy2*cosd(lineangle);
        c1new.accel(1)=c1new.accel(1)-(ax1*cosd(lineangle)+ay1*sind(lineangle));
        c1new.accel(2)=c1new.accel(2)-(-ax1*sind(lineangle)+ay1*cosd(lineangle));
        c2new.accel(1)=c2new.accel(1)+(ax1*cosd(lineangle)+ay1*sind(lineangle));
        c2new.accel(2)=c2new.accel(2)+(-ax1*sind(lineangle)+ay1*cosd(lineangle));
    catch
        disp('in the circles function')
        pause;
    end

    
endfunction
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
polygon.points=[
73.6 0
95.8 0
72.8 -9.3
94.4 -15.2];
polygon.lines=[];
//polygon.lines=[polygon.lines; polygon.points(1,:) polygon.points(2,:)];
polygon.lines=[polygon.lines; polygon.points(1,:) polygon.points(3,:)];
polygon.lines=[polygon.lines; polygon.points(2,:) polygon.points(4,:)];//setup an open bin
polygon.lines=[polygon.lines; polygon.points(4,:) polygon.points(3,:)];

circlecount=5;

c.mass=10;
c.radius=1;
c.pos=[mean(polygon.points(:,1)) mean(polygon.points(:,2))];
c.vel=[0 0]';
c.accel=[0 0];



circle(1:circlecount)=c;
for i=1:size(circle,1)
    circle(i).pos=[circle(i).pos(1)+2*c.radius*grand(1, 1, "unf", -2, 2) circle(i).pos(2)+3*c.radius*i];
end

dt=0.03;
Time1=10;
scf();        
elasticity=0.1;

draw=1;
disp('Running...')
for t=0:dt:Time1
    if draw==1
        drawlater();
        clf;
        
        for i=1:size(polygon.lines,1)
            plot(polygon.lines(i,[1,3]),polygon.lines(i,[2,4]),'-');
        end
        
    end
    stopped=0;
    for i=1:size(circle,1)
        if  ((abs(circle(i).vel(2)))<0.5 & t > 1) //circle(i).pos(2)<c.radius*0.5  
            stopped=stopped+1;
        end
        circle(i)=gravityforce(circle(i));
        circle(i)=dampingforce(circle(i));
        circle_old=circle(i);
        circle(i)=new_pos_vel(circle(i),dt);
        circle(i)=findCircDistances(circle(i));
        if  circle(i).inside_polygon==0
            //pause
            indexofline=find(circle(i).linedist<circle(i).radius);
            indexofline=min(indexofline);
            linetheta=atand((polygon.lines(indexofline,4)-polygon.lines(indexofline,2))/(polygon.lines(indexofline,3)-polygon.lines(indexofline,1)));
            circle(i)=getnewspeed(circle(i),circle_old,linetheta,elasticity);
            //circle(i).pos(2)=circle_old.pos;
        elseif circle(i).clearofcircles==0
            
            indexofcircle=find(circle(i).circdist(find(circle(i).circdist>0))<2*circle(i).radius);
            if indexofcircle>=i
                indexofcircle=indexofcircle+1;
            end
            indexofcircle=min(indexofcircle);
            [circle(i) circle(indexofcircle)]=getnewspeedcircles(circle(i),circle_old,circle(indexofcircle),elasticity);
            //pause
            //circle(i).pos=circle_old.pos;
        end
        if draw==1
            disp([t i circle(i).pos circle(i).vel' circle(i).accel circle(i).clearofcircles circle(i).inside_polygon])
            draw_circle(circle(i));
        end
    end
    if stopped==circlecount
        abort
    end
    
    if draw==1
        h=gca();h.axes_visible=['off','off','off'];h.box='off';//p=max(h.data_bounds);q=min(h.data_bounds);h.data_bounds=[q q;p p];
        drawnow();
    end
end
 if draw==0 
    
    for i=1:size(polygon.lines,1)
        plot(polygon.lines(i,[1,3]),polygon.lines(i,[2,4]),'-');
    end
    for i=1:size(circle,1)
        draw_circle(circle(i));
    end
    h=gca();h.axes_visible=['off','off','off'];h.box='off';p=max(h.data_bounds);q=min(h.data_bounds);h.data_bounds=[q q;p p];
end
disp('completed.')
