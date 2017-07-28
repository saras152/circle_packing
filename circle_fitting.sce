 
clear;
clearglobal
function draw_circle(c)
    theta=linspace(0,2*%pi,100);
    x=c.pos(1)+c.radius*cos(theta);
    y=c.pos(2)+c.radius*sin(theta);
    xstring(c.pos(1),c.pos(2),string(i));
    plot(x,y);
endfunction
//function [c]=circlerepel(c,c_source)
////    //disp('in the circle function.')
//////pause;
////    Constant=1e1;
////    c1=complex(c.pos(1),c.pos(2));
////    for i=1:size(c_source,1)
////        c2=complex(c_source(i).pos(1),c_source(i).pos(2));
////        c2=c1-c2;
////        dist=abs(c2);
////        dire=atan(imag(c2),real(c2));
////        if  dist==0
////            accel=0;
////        elseif dist <= (c.radius+c_source(i).radius)
////            if dist<(c.radius+c_source(i).radius)
////            accel=Constant*1*c_source(i).mass*exp(%i*dire)/(c.radius+c_source(i).radius)^2;
////        else
////            accel=Constant*c_source(i).mass*exp(%i*dire)/(dist)^2;
////        end
////        c.accel=c.accel+[real(accel) imag(accel)];
////    end
//endfunction
//
//function c=linerepel(c,polyg)
////disp('in the lines function.')
////pause;
//    segments=10;
//    Constant=1e1;
//    c1=complex(c.pos(1),c.pos(2));
//    corners=polyg.points;
//    linearray=polyg.lines;
//    corners($+1,:)=corners(1,:);
//    angle1=0;
//    dist=[];
//    for i=1:(size(corners,1)-1)
//        c2=complex(corners(i,1),corners(i,2)) - c1;//reference to the point
//        c3=complex(corners(i+1,1),corners(i+1,2)) - c1;//reference to the point
//        c23=c3-c2;//line that connects the two points
//        a2=(atand(imag(c2),real(c2)));
//        a3=(atand(imag(c3),real(c3)));
//
//        angle1=angle1+a3-a2;
//        if a3-a2>180
//            angle1=angle1-360;
//        elseif +a3-a2< -180
//            angle1=angle1+360;
//        end
//        x1=corners(i,1);y1=corners(i,2);
//        x2=corners(i+1,1);y2=corners(i+1,2);
//        x0=c.pos(1);y0=c.pos(2);
//        d=abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2);
//        dist=[dist; d];
//        //disp(dist)
//    end
//    if angle1<0.01 | min(dist) <= c.radius
//        inside_polygon=0;
//    else
//        inside_polygon=1;
//    end
//    c.inside_polygon=inside_polygon;
//    for i=1:size(linearray,1)
//        if inside_polygon==0
//            if dist==c.radius//single point of contact
//                
//            else//less = two points of contact. 
//                
//            end
//            
//        end        
//                
//    end
//    
////    if  inside_polygon == 1
////        for i=1:size(linearray,1)
////            lenline=sqrt((linearray(i,3)-linearray(i,1))^2+(linearray(i,4)-linearray(i,2))^2);
////            dl=lenline/segments;
////            for j=1:segments
////                xmid=(linearray(i,1)*(j-0.5) + linearray(i,3)*(segments-j+0.5))/segments;
////                ymid=(linearray(i,2)*(j-0.5) + linearray(i,4)*(segments-j+0.5))/segments;
////                c2=complex(xmid,ymid);
////                c2=c1-c2;
////                dist=abs(c2);
////                dire=atan(imag(c2),real(c2));
////                if  dist==0
////                    accel=0;
////                elseif dist < (c.radius)
////                    accel=Constant*1000*dl*exp(%i*dire)/(c.radius)^2;
////                else
////                    accel=Constant*dl*exp(%i*dire)/(dist)^2;
////                end
////                c.accel=c.accel+[real(accel) imag(accel)]
////            end
////        end
////    else
////        accel=0;
////    end
//endfunction

function c=dampingforce(c)
    k=0.9;pow=1;
    c.accel=c.accel-k*(c.vel^pow);
endfunction
function c=gravityforce(c)
    gravity.accel=[0 -9.8];
    c.accel=gravity.accel;
endfunction
function c=update_pos_vel(c,dt)
    
    global polygon
    global circle
    c1.vel=c.vel+c.accel*dt;
    c1.pos=c.pos+c.vel*dt;
    //pause
    cnew=complex(c1.pos(1),c1.pos(2));
    corners=polygon.points;
    linearray=polygon.lines;
    corners($+1,:)=corners(1,:);
    angle1=0;
    dist=[];
    for i=1:(size(corners,1)-1)
        c2=complex(corners(i,1),corners(i,2)) - cnew;//reference to the point
        c3=complex(corners(i+1,1),corners(i+1,2)) - cnew;//reference to the point
        c23=c3-c2;//line that connects the two points
        a2=(atand(imag(c2),real(c2)));
        a3=(atand(imag(c3),real(c3)));

        angle1=angle1+a3-a2;
        if a3-a2 > 180
            angle1=angle1-360;
        elseif a3-a2 < -180
            angle1=angle1+360;
        end
        x1=corners(i,1);y1=corners(i,2);
        x2=corners(i+1,1);y2=corners(i+1,2);
        x0=c1.pos(1);y0=c1.pos(2);
        d=abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2);
        dist=[dist; d];
        //disp(dist)
    end
    if angle1<0.01 | min(dist) <= c.radius
        inside_polygon=0;
    else
        inside_polygon=1;
    end
    distc=100*ones(size(circle,1),size(circle,1));
    
    for i=1:size(circle,1)
        for j=i+1:size(circle,1)
            distc(i,j)=[sqrt((circle(i).pos(1)-circle(j).pos(1))^2 + (circle(i).pos(2)-circle(j).pos(2))^2)];
        end
    end
    if  min(distc) <= 2*c.radius
        overlapping_circles=1;
    else
        overlapping_circles=0;
    end
    //disp(dist)
    
    elasticity=0.5;
    if  inside_polygon==0
        c.vel=0;
        lin1=find(dist<=c.radius);
        
    elseif overlapping_circles==1
        //index=
        c.vel=0;
        
    else
        c.vel=c1.vel;
        c.pos=c1.pos;
    end
endfunction
global polygon
global circle
polygon.points=[
0 0
30 0
30 30
0 30];
polygon.lines=[];
polygon.lines=[polygon.lines; polygon.points(1,:) polygon.points(2,:)];
polygon.lines=[polygon.lines; polygon.points(2,:) polygon.points(3,:)];
polygon.lines=[polygon.lines; polygon.points(3,:) polygon.points(4,:)];//setup an open bin
polygon.lines=[polygon.lines; polygon.points(4,:) polygon.points(1,:)];

circlecount=2;

c.mass=10;
c.radius=1;
c.pos=[mean(polygon.points(:,1)) mean(polygon.points(:,2))];
c.vel=[0 0];
c.accel=[0 0];
c.inside_polygon=1;


circle(1:circlecount)=c;
for i=1:size(circle,1)
    circle(i).pos=[circle(i).pos(1)+2*rand(1) circle(i).pos(2)+2*c.radius*grand(1, 1, "uin", 1, 10)];
end



dt=0.01;
Time1=10;
scf(100);
for t=0:dt:Time1
    
    drawlater();
    clf;
    
    for i=1:size(polygon.lines,1)
        plot(polygon.lines(i,[1,3]),polygon.lines(i,[2,4]),'-');
    end
    for i=1:size(circle,1)
        //disp('in the main function.')

        if circle(i).pos(2)<c.radius*0.5 
            abort;
        end
        //pause
        circle(i)=gravityforce(circle(i));
        //circle(i)=circlerepel(circle(i),circle);
        //circle(i)=linerepel(circle(i),polygon);
        circle(i)=dampingforce(circle(i));
        circle(i)=update_pos_vel(circle(i),dt);
        disp([t i circle(i).pos circle(i).vel circle(i).accel ])
        //plot(circle(i).pos(1),circle(i).pos(2),'*');
        draw_circle(circle(i));
        h=gca();h.axes_visible=['off','on','off'];h.box='off';p=max(h.data_bounds);q=min(h.data_bounds);h.data_bounds=[q q;p p];
        //pause;
    end
    drawnow();
end


