clear
function draw_circle(c)
    theta=linspace(0,2*%pi,16);
    x=c.pos(1)+c.radius*cos(theta);
    y=c.pos(2)+c.radius*sin(theta);
    xstring(c.pos(1),c.pos(2),string(i));
    plot(x,y);
endfunction
polygon.points=[
73.6 0
95.8 0
94.4 -15.2
72.8 -9.3
];
polygon.lines=[];
polygon.lines=[polygon.lines; polygon.points(1,:) polygon.points(2,:)];
polygon.lines=[polygon.lines; polygon.points(2,:) polygon.points(3,:)];
polygon.lines=[polygon.lines; polygon.points(3,:) polygon.points(4,:)];
polygon.lines=[polygon.lines; polygon.points(4,:) polygon.points(1,:)];

circlecount=5;
c.radius=1;
polygon.lines=[polygon.lines;polygon.lines(1,:)];
for i=1:size(polygon.lines,1)-1
    linetheta1=(atand((polygon.lines(i,4)-polygon.lines(i,2))/(polygon.lines(i,3)-polygon.lines(i,1))));
    linetheta2=(atand((polygon.lines(i+1,4)-polygon.lines(i+1,2))/(polygon.lines(i+1,3)-polygon.lines(i+1,1))));
    d=(c.radius/sind((linetheta2-linetheta1)/2));
    circle(i).pos(1)=polygon.lines(i,3)+d*cosd((linetheta2+linetheta1)/2);
    circle(i).pos(2)=polygon.lines(i,4)+d*sind((linetheta2+linetheta1)/2);
    circle(i).radius=c.radius;
end
scf;
for i=1:size(polygon.lines,1)-1
    plot(polygon.lines(i,[1,3]),polygon.lines(i,[2,4]),'-');
end
for i=1:size(circle,1)
    draw_circle(circle(i));
end  
