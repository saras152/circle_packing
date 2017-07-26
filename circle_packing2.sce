// exec('C:\Users\raghunath\Desktop\circle_packing2.sce',-1)
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
function [cn,on] = inpoly(p,node,myedge,TOL)
    //    (c) 2006 Darren Engwirda//    
    //    // Geometry
    //dtheta = %pi/15;
    //theta  = (-%pi:dtheta:(%pi-dtheta))';
    //node   = [cos(theta) sin(theta)];
    //n      = size(node,1);
    //cnect  = [(1:n-1)' (2:n)'; n 1];
    //p = 3*(rand(10000,2)-0.5);
    //p1=node;
    //p2=[
    //-0.5 -0.5
    //0.5 -0.5
    //0.5 0.5
    //-0.5 0.5 
    //];
    //n1=size(p1,1);
    //n2=size(p2,1);
    //c1 = [(1:n1-1)', (2:n1)'; n1, 1];
    //c2 = [(1:n2-1)', (2:n2)'; n2, 1];
    //node=[p1;p2];
    //cnect=[
    //c1
    //c2+n1];
    //tic, in = inpoly(p,node,cnect); t1 = toc;
    //scf;plot(p(find(in==1),1),p(find(in==1),2),'b.');plot(p(find(in==0),1),p(find(in==0),2),'r.');

    nargin=argn(2);
    if nargin<4
        TOL = 1.0e-12;
        if nargin<3
            myedge = [];
            if nargin<2
                error('Insufficient inputs');
            end
        end
    end
    nnode = size(node,1);
    if isempty(myedge)                                                           // Build edge if not passed
        myedge = [(1:nnode-1)' (2:nnode)'; nnode 1];
    end
    if size(p,2)~=2
        error('P must be an Nx2 array.');
    end
    if size(node,2)~=2
        error('NODE must be an Mx2 array.');
    end
    if size(myedge,2)~=2
        error('EDGE must be an Mx2 array.');
    end
    if max(myedge(:))>nnode | or(myedge(:)<1)
        error('Invalid EDGE.');
    end

    //// PRE-PROCESSING
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    n  = size(p,1);
    nc = size(myedge,1);

    // Choose the direction with the biggest range as the "y-coordinate" for the
    // test. This should ensure that the sorting is done along the best
    // direction for long and skinny problems wrt either the x or y axes.
    dxy = max(p,'r')-min(p,'r');
    if dxy(1)>dxy(2)
        // Flip co-ords if x range is bigger
        p = p(:,[2,1]);
        node = node(:,[2,1]);
    end
    tol = TOL*min(dxy);

    // Sort test points by y-value
    [y,i] = gsort(p(:,2),'g','i');
    x = p(i,1);

    //// MAIN LOOP
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cn = zeros(n,1);     // Because we're dealing with mod(cn,2) we don't have
    // to actually increment the crossing number, we can
    // just flip a logical at each intersection (faster!)
    on = cn;
    for k = 1:nc         // Loop through edges

        // Nodes in current edge
        n1 = myedge(k,1);
        n2 = myedge(k,2);

        // Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
        //           - also get xmin = min(x1,x2), xmax = max(x1,x2)
        y1 = node(n1,2);
        y2 = node(n2,2);
        if y1<y2
            x1 = node(n1,1);
            x2 = node(n2,1);
        else
            yt = y1;
            y1 = y2;
            y2 = yt;
            x1 = node(n2,1);
            x2 = node(n1,1);
        end
        if x1>x2
            xmin = x2;
            xmax = x1;
        else
            xmin = x1;
            xmax = x2;
        end

        // Binary search to find first point with y<=y1 for current edge
        if y(1)>=y1
            start = 1;
        elseif y(n)<y1
            start = n+1;       
        else
            lower = 1;
            upper = n;
            for j = 1:n
                start = round(0.5*(lower+upper));
                if y(start)<y1
                    lower = start;
                elseif y(start-1)<y1
                    break;
                else
                    upper = start;
                end
            end
        end
        // Loop through points
        for j = start:n
            // Check the bounding-box for the edge before doing the intersection
            // test. Take shortcuts wherever possible!

            Y = y(j);   // Do the array look-up once & make a temp scalar
            if Y<=y2
                X = x(j);   // Do the array look-up once & make a temp scalar
                if X>=xmin
                    if X<=xmax
                        //pause;
                        // Check if we're "on" the edge
                        on(j) = on(j) | (abs((y2-Y)*(x1-X)-(y1-Y)*(x2-X))<tol);
                        // Do the actual intersection test
                        if (Y<y2) & ((y2-y1)*(X-x1)<(Y-y1)*(x2-x1))
                            cn(j) = ~cn(j);
                        end

                    end
                elseif Y<y2   // Deal with points exactly at vertices
                    // Has to cross edge
                    cn(j) = ~cn(j);
                end
            else
                // Due to the sorting, no points with >y
                // value need to be checked
                break
            end
        end

    end

    // Re-index to undo the sorting
    cn(i) = cn|on;
    on(i) = on;

endfunction      // inpoly()

function circlenodes=movecircle(circlenodes,CircleID,dx,dy,NodesPerCircle)
    //    pause

    noderange=1:NodesPerCircle;
    noderange=(CircleID-1)*NodesPerCircle + noderange;
    circlenodes(noderange,1)=circlenodes(noderange,1)+dx;
    circlenodes(noderange,2)=circlenodes(noderange,2)+dy;

endfunction

function out=mysign(val)
    if val>=0
        out=1;
    else
        out=-1;
    end
endfunction

function drawpolygons(polynodes,polyconnections)
    newsequence=[];
    for i=1:size(polyconnections,1)
        newsequence=[newsequence;polyconnections(i,1);polyconnections(i,2);0];
    end
    newpoints=[];
    //pause
    for i=1:length(newsequence)
        if newsequence(i)==0
            newpoints=[newpoints; %nan %nan];
        else
            newpoints=[newpoints;polynodes(newsequence(i),:)];

        end
    end
    plot(newpoints(:,1),newpoints(:,2),'r-');
endfunction
function dist=boundarydistance(point,boundary)
    boundary=[boundary;boundary(1,:)];
    dist=orthProj(boundary,point);
endfunction
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
draw=1;
processfurther=1;
if  processfurther==1
    tic;
    polynodes=[
    84.991668 3.75
    //112.157 3.75 
    113.656 5.32
    113.378 9.56
    111.687 10.922
    84.754 7.3757
    ];
    np      = size(polynodes,1);
    polyconnections  = [(1:np-1)' (2:np)'; np 1];
    onlypolynodes=polynodes;
    circlecount=200;
    radius=0.8128/2;
    r=radius*(1+1e-4);
    NodesPerCircle=16;
    theta=(1:NodesPerCircle)'*2*%pi/NodesPerCircle;

    circlenodeinit=np;
    circlenodes=[];
    circleconnections=[];

    //bottom left point
    x=polynodes(:,1);
    y=polynodes(:,2);
    x=x-mean(x);y=y-mean(y);
    indices_q1=intersect(find(x>0), find(y>0));
    indices_q3=intersect(find(x<0), find(y<0));
    indices_q=indices_q3;
    dist_cg=sqrt(x.^2+y.^2);
    farthest_index=indices_q(min(find(dist_cg(indices_q)==max(dist_cg(indices_q)))));
    farthest_pt=polynodes(farthest_index,:);
    lineindex1=find(polyconnections(:,1)==farthest_index);
    lineindex2=find(polyconnections(:,2)==farthest_index);
    linetheta1=((atand((polynodes(polyconnections(lineindex1,2),2)-polynodes(polyconnections(lineindex1,1),2)),(polynodes(polyconnections(lineindex1,2),1)-polynodes(polyconnections(lineindex1,1),1)))));
    linetheta2=((atand((polynodes(polyconnections(lineindex2,1),2)-polynodes(polyconnections(lineindex2,2),2)),(polynodes(polyconnections(lineindex2,1),1)-polynodes(polyconnections(lineindex2,2),1)))));
    d=abs(r/sind((linetheta2-linetheta1)/2));
    pos=farthest_pt  + [mysign(farthest_pt(1))*d*cosd((linetheta2+linetheta1)/2) mysign(farthest_pt(2))*d*sind((linetheta2+linetheta1)/2) ];
    //    seedpos(1)=pos(1)-2*r*cosd(linetheta1);
    //    seedpos(2)=pos(2)-2*r*sind(linetheta1);
    circleCoordinates=[];
    oldpos=[];
    forward=1;
    if draw==1
        clf;
        drawlater();
        drawpolygons(polynodes,polyconnections);
        drawnow();
    end
    dseedtheta=1;
    radiusfactor=1;
    narrowpassage=0;
    //disp('##starting here. ')
    for i=1:circlecount
        inside=0;
        seedtheta=0;
        radiusfactor=1;
        if i==137
            //            pause
        end
        //        pause
        //disp('planning to place the wire number : ' + string(i))
        while(inside==0)
            thiscirclenode=[pos(1)+radius*cos(theta) pos(2)+radius*sin(theta)];
            nc= size(thiscirclenode,1);
            thiscircleconnections=circlenodeinit + [(1:nc-1)' (2:nc)'; nc 1];
            [in,on]=inpoly(thiscirclenode,polynodes,polyconnections)

            if length(find(in==0))==0 & length(find(on==1))==0//perfectly inside.
                inside=1;
                //pos=pos+r*0.01*rand(1,2);
                circleCoordinates=[circleCoordinates;pos]
                oldpos=pos;
                //disp('found the fitting place for the wire ' + string(i))
                if seedtheta<>0
                    //disp('seedtheta is non zero. flipping the direction')
                    seedtheta=0;
                    if boundarydistance(pos,onlypolynodes) < 1.5*r
                        if narrowpassage==0
                            forward=-forward;
                        else
                            narrowpassage=0;
                        end
                        
                    end
                end
                pos=pos + [ forward*2*r*cosd(linetheta1+seedtheta)  forward*2*r*sind(linetheta1+seedtheta)] //+r*0.001*grand(1,2,'unf',-1,1)
//            elseif length(find(in==0))==0 & length(find(on==0))==0//overlapping. may not happen for now.
//                searching=1;
//                forward=-forward;
                //disp('found overlapping position.')
                //pos=pos + [ 2*r*cosd(linetheta1)  2*r*sind(linetheta1)]
            else
                if length(find(in==0))==0 & length(find(on==0))==0//overlapping
                    disp('narrowpassage activated')
                    narrowpassage=1;
                    forward=-forward;
                        radiusfactor=radiusfactor+1;
                    seedtheta=0+1*forward*dseedtheta/radiusfactor
                end
                
                if length(find(in==1))==0//perfectly outside
                    //                    if radiusfactor>1
                    //                        seedtheta=seedtheta-forward*dseedtheta;
                    //                    end
                    if narrowpassage==1
                        radiusfactor=radiusfactor+1;
                        seedtheta=0+1*forward*dseedtheta/radiusfactor;
                    else
                         
                        seedtheta=seedtheta+1*forward*dseedtheta/radiusfactor;
                        //disp('FF')
                    end
//                    pause
                else
                    seedtheta=seedtheta+forward*dseedtheta/radiusfactor;
                end


                if abs(seedtheta) > 180
                    //                    pause;
                    radiusfactor=radiusfactor+1;
                    seedtheta=seedtheta-forward*180/radiusfactor;
                    disp('increasing radius factor')
                    //                    if boundarydistance(pos,onlypolynodes) < 1.5*r
                    //                        forward=-forward;
                    //                    end
                end
                pos=oldpos + [ forward*radiusfactor*2*r*cosd(linetheta1+seedtheta) forward*radiusfactor*2*r*sind(linetheta1+seedtheta)] //+r*0.05*rand(1,2);
            end
            if draw==1
                drawlater();
                clf;
                drawpolygons([polynodes;thiscirclenode],[polyconnections;thiscircleconnections]);
                drawnow();
            end
        end
        polynodes=[polynodes;thiscirclenode];
        polyconnections=[polyconnections;thiscircleconnections];
        circlenodeinit=circlenodeinit+size(thiscirclenode,1);
        if draw==1
            drawlater();
            clf;
            drawpolygons(polynodes,polyconnections);
            drawnow();
        end






        //        circlenodes=[circlenodes;thiscirclenode];
        //        circleconnections=[circleconnections;thiscircleconnections];
    end
    //filledcircles=0;
    ////for i=1:circlecount
    ////    circlenodes=movecircle(circlenodes,i,20*rand(1,1),10*rand(1,1),16);
    ////end
    ////totalnodes=[polynodes;circlenodes];
    ////totalconnections=[polyconnections;circleconnections];
    ////p = 120*(rand(100000,2)+0);
    ////in = inpoly(p,totalnodes,totalconnections); 
    ////scf;
    ////plot(p(find(in==1),1),p(find(in==1),2),'b.');
    ////plot(p(find(in==0),1),p(find(in==0),2),'r.');
    //////plot(totalnodes(:,1),totalnodes(:,2),'g')
    //
    //while(filledcircles<circlecount)
    //    thiscircle=circlenodes(1:16);
    //    thiscircleconnections=circleconnections(1,16);
    //end
end
execstr('tm=toc();disp(''Simulation of the circle packing took ''+string(round(tm))+'' seconds.'');')
