function plot_cube(origin,length,color)
% Define the vertexes of the unit cubic
ver = [1 1 0;
    0 1 0;
    0 1 1;
    1 1 1;
    0 0 1;
    1 0 1;
    1 0 0;
    0 0 0];

%  Define the faces of the unit cubic
fac = [1 2 3 4;
    4 3 5 6;
    6 7 8 5;
    1 2 8 7;
    6 7 1 4;
    2 3 5 8];
cube = [ver(:,1)*length(1),ver(:,2)*length(2),ver(:,3)*length(3)];
cube = [cube(:,1)+origin(1),cube(:,2)+origin(2),cube(:,3)+origin(3)];
patch('Faces',fac,'Vertices',cube,'FaceVertexCData',color,'FaceColor','flat');
end