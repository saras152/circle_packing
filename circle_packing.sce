//          exec('C:\Users\raghunath\Desktop\circle_packing.sce',-1)
clear
clearglobal
tic
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////FUNCTION DEFINITIONS START FROM HERE//////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function inside = point_in_polygon(c)

    global polygon;
    global radius;
    xpol=polygon.points(:,1)';
    ypol=polygon.points(:,2)';
    xpoint=c.pos(1);
    ypoint=c.pos(2);
    npol = size(xpol,'*')
    inside = 0
    j = npol; // j is the previous vertice 
    i = 1
    while  i <= npol
        if ((((ypol(i) <= ypoint) & (ypoint < ypol(j)))|((ypol(j) <= ypoint) & (ypoint < ypol(i)))) & ..
           (xpoint < ((xpol(j) - xpol(i))/(ypol(j) - ypol(i))) * (ypoint - ypol(i)) + xpol(i)))
              inside = 1-inside;
        end
        i = i + 1;
        j = i - 1;
    end
endfunction
function c=findCircDistances(cold)
    global polygon
    global circle
    c=cold;
//    linedist=[];
//    circdist=zeros(size(circle,1),1);
//    x0=c.pos(1);y0=c.pos(2);
//    for j=1:size(polygon.lines,1)
//        x1=polygon.lines(j,1);y1=polygon.lines(j,2);
//        x2=polygon.lines(j,3);y2=polygon.lines(j,4);
//        linedist=[linedist abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2)];
//    end
//    c.linedist=linedist;
    pgn=[polygon.points;polygon.points(1,:)];
    //disp(orthProj(pgn,c.pos))
    if   point_in_polygon(c)==0  | orthProj(pgn,c.pos) <= radius//| min(linedist) < radius 
        c.inside_polygon=0;
        //disp('outside polygon')
        //pause
    else
        c.inside_polygon=1;
    end
    if size(circle,1) >1
    
        for j=1:size(circle,1)
            circdist(j)=[sqrt((c.pos(1)-circle(j).pos(1))^2 + (c.pos(2)-circle(j).pos(2))^2)];
        end
        c.circdist=circdist;
        if  min(circdist(find(circdist>0))) < 2*radius
            c.clearofcircles=0; 
            //disp('overlapping circles')
            //pause
        else
            c.clearofcircles=1;
        end
    else
        c.clearofcircles=1;
    end
    
endfunction
function out=mysign(val)
    if val>=0
        out=1;
    else
        out=-1;
    end
endfunction
function draw_circle(c)
    global radius;
    theta=linspace(0,2*%pi,16);
    x=c.pos(1)+radius*cos(theta);
    y=c.pos(2)+radius*sin(theta);
    xstring(c.pos(1),c.pos(2),string(i));
    plot(x,y);
endfunction
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////ACTUAL CODE STARTS FROM HERE//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

polygon.points=[
84.991668 3.75
112.157 3.75 
113.656 5.32
113.378 9.56
111.687 10.922
84.754 7.3757
];
polygon.lines=[];
for i=1:size(polygon.points,1)-1
    polygon.lines=[polygon.lines; polygon.points(i,:) polygon.points(i+1,:)];
end
polygon.lines=[polygon.lines; polygon.points(i+1,:) polygon.points(1,:)];

//polygon.lines=[polygon.lines; polygon.points(1,:) polygon.points(2,:)];
//polygon.lines=[polygon.lines; polygon.points(2,:) polygon.points(3,:)];
//polygon.lines=[polygon.lines; polygon.points(3,:) polygon.points(4,:)];
//polygon.lines=[polygon.lines; polygon.points(4,:) polygon.points(1,:)];
circlecount=11*12;
radius=0.8128/2;
//Right Most Point
x=polygon.points(:,1);
y=polygon.points(:,2);
x=x-mean(x);y=y-mean(y);
indices_q1=intersect(find(x>0), find(y>0));
indices_q=indices_q1;
dist_cg=sqrt(x.^2+y.^2);
farthest_index=indices_q(min(find(dist_cg(indices_q)==max(dist_cg(indices_q)))))
farthest_pt=polygon.points(farthest_index,:);
lineindex1=vectorfind(polygon.lines(:,1:2),farthest_pt);
lineindex2=vectorfind(polygon.lines(:,3:4),farthest_pt);
linetheta1=abs((atand((polygon.lines(lineindex1,4)-polygon.lines(lineindex1,2))/(polygon.lines(lineindex1,3)-polygon.lines(lineindex1,1)))));
linetheta2=abs((atand((polygon.lines(lineindex2,4)-polygon.lines(lineindex2,2))/(polygon.lines(lineindex2,3)-polygon.lines(lineindex2,1)))));
d=abs(radius/sind((linetheta2-linetheta1)/2));
pos=farthest_pt  - [mysign(farthest_pt(1))*d*cosd((linetheta2+linetheta1)/2) mysign(farthest_pt(2))*d*sind((linetheta2+linetheta1)/2) ];
steps=20;
dx=(max(polygon.points(:,1))-min(polygon.points(:,1)))/steps;
dy=(max(polygon.points(:,2))-min(polygon.points(:,2)))/steps;
circle=[];
global circle;
global polygon;
global radius;
scf;
draw=1;
tic
for i=1:circlecount
    circle(i).pos=pos;  
    circle(i).mobilex=0;
    circle(i).mobiley=1;
//    pause;
    iterationsmax=15;
    iterationcounter=1;
    tic()
    while circle(i).mobilex==1 | circle(i).mobiley==1
        cold=circle(i);
        fraction=iterationsmax/(iterationcounter)^.9;
        //fraction=iterationsmax-iterationcounter+1
        if circle(i).mobilex==1
            if iterationcounter < 0.4*iterationsmax | iterationcounter >= 0.5*iterationsmax //| iterationcounter < 0.4*iterationsmax | iterationcounter >= 0.6*iterationsmax
                circle(i).pos(1)=circle(i).pos(1) - mysign(farthest_pt(1))*dx*0.5*(fraction)*grand(1,1,'unf',-0.1,1);
            else
                circle(i).pos(1)=circle(i).pos(1) - mysign(farthest_pt(1))*dx*2*(fraction)*grand(1,1,'unf',-0.1,1);
            end
        end
        if circle(i).mobiley==1
            if iterationcounter < 0.4*iterationsmax | iterationcounter >= 0.5*iterationsmax //| iterationcounter < 0.4*iterationsmax | iterationcounter >= 0.6*iterationsmax
                circle(i).pos(2)=circle(i).pos(2) - mysign(farthest_pt(2))*dy*0.5*(fraction)*grand(1,1,'unf',-0.1,1);
            else
                circle(i).pos(2)=circle(i).pos(2) - mysign(farthest_pt(2))*dy*2*(fraction)*grand(1,1,'unf',-0.1,1);//
            end
            
        end
        circle(i)=findCircDistances(circle(i));
        
        if circle(i).clearofcircles==0 | circle(i).inside_polygon==0
            if circle(i).mobilex==1
                circle(i)=cold;
                circle(i).mobilex=0;
                circle(i).mobiley=1;
                //disp('x made to 0 for ' + string(i))
                //pause
            elseif circle(i).mobiley==1
                circle(i)=cold;
                circle(i).mobiley=0;
                circle(i).mobilex=1;
                //disp('y made to 0 for ' + string(i))
               // pause
                
                iterationcounter=iterationcounter+1;
            end
            if iterationcounter>iterationsmax
                circle(i)=cold;
                circle(i).mobiley=0;
                circle(i).mobilex=0;
                //disp('both stopped for ' + string(i))
                //pause
                
            end
        end
    end
    
    execstr('tm=toc();disp(''finding the location took ''+string((tm))+'' seconds.'');')
    if draw==1
        drawlater();
        clf;
        
        for i=1:size(polygon.lines,1)
            plot(polygon.lines(i,[1,3]),polygon.lines(i,[2,4]),'-');
        end
        for i=1:size(circle,1)
            draw_circle(circle(i));
        end  
        drawnow();
    end
end

if draw==0
    drawlater();
    clf;
    
    for i=1:size(polygon.lines,1)
        plot(polygon.lines(i,[1,3]),polygon.lines(i,[2,4]),'-');
    end
    for i=1:size(circle,1)
        draw_circle(circle(i));
    end  
    drawnow();
end

 //execstr('tm=toc();disp(''Simulation of the circle packing took ''+string(round(tm))+'' seconds.'');tic();')


