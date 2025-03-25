function [v_grid, f_grid, c_grid, i_grid, x_grid] = generate_triangular_grid(rows, cols, l, ox, oy)
% Generate the triangular grids given the number of rows and cols, and the
% length of a unit, offset distance of x and y.

% Default parameter values
if nargin < 3, l = 1; end
if nargin < 4, ox = 0; end
if nargin < 5, oy = 0; end

% Initialize variables
v_grid = []; % vertices
f_grid = []; % faces:connections of vertices
c_grid = []; % centroid of faces
i_grid = []; % Orietation of faces, 0:upwards, 1:downwards
x_grid = {}; % Neighboring triangle IDs

h = sqrt(3)/2 * l; % Height of triangle
idE = 0;

% Generate vertices
for row = 0:rows
    for col = 0:cols
        if mod(row,2)==0
            v_grid(end+1,:) = [l*col+ox, h*row+oy, 0];
        else
            v_grid(end+1,:) = [l*col+l/2+ox, h*row+oy, 0];
        end
    end
end

% Generate faces and centroids
for row = 0:rows-1
    iList = [];
    for col = 0:cols-1
        if mod(row,2)==0 % odd row
            m=1; n=0;
            iList=[iList,0,1];
        else % even row
            m=0; n=1;
            iList=[iList,1,0];
        end

        % Triangle upwards
        % Vertices indice for triangles
        n0=col+row*(cols+1);
        n1=col+row*(cols+1)+1;
        n2=col+(row+1)*(cols+1)+n;

        fL1=[n0+1,n1+1,n2+1];

        cL1=[mean(v_grid(fL1(1:2),1)), v_grid(fL1(1),2)+sqrt(3)/6*l, 0]; % centroid upward

        % Triangle downwards
        % Vertices indice for triangles
        n0=col+row*(cols+1)+m;
        n1=col+(row+1)*(cols+1);
        n2=col+(row+1)*(cols+1)+1;


        fL2=[n0+1,n1+1,n2+1];

        cL2=[mean(v_grid(fL2(2:3),1)), v_grid(fL2(2),2)-sqrt(3)/6*l, 0]; % centroid downward

        if mod(row,2)==0
            f_grid=[f_grid;fL1;fL2];
            c_grid=[c_grid;cL1;cL2];

            % Neighbor IDs (odd row)
            if row==0 % first row
                if col==0 % first col
                    x_grid{end+1}=idE+2;
                else
                    x_grid{end+1}=[idE,idE+2];
                end
            else 
                if col==0
                    x_grid{end+1}=[idE+2,idE-(cols-1)*2+1];
                else
                    x_grid{end+1}=[idE,idE-(cols-1)*2+1,idE+2];
                end
            end

            idE=idE+1;

            if row==rows-1 % final row
                if col==cols-1
                    x_grid{end+1}=idE;
                else
                    x_grid{end+1}=[idE,idE+2];
                end
            else
                if col==cols-1
                    x_grid{end+1}=[idE,idE+(cols-1)*2];
                else
                    x_grid{end+1}=[idE,idE+(cols-1)*2,idE+2];
                end
            end

            idE=idE+1;

        else % even rows
            f_grid=[f_grid;fL2;fL1]; % reversed order triangles
            c_grid=[c_grid;cL2;cL1];

            if row==rows-1 % final row
                if col==0
                    x_grid{end+1}=idE+2;
                else
                    x_grid{end+1}=[idE,idE+2];
                end
            else
                if col==0
                    x_grid{end+1}=[idE+2,idE+(cols-1)*2];
                else
                    x_grid{end+1}=[idE,idE+(cols-1)*2,idE+2];
                end
            end

            idE=idE+1;

            if row==0 % first row
                if col==cols-2
                    x_grid{end+1}=idE;
                else
                    x_grid{end+1}=[idE,idE+2];
                end
            else
                if col==cols-2
                    x_grid{end+1}=[idE,idE-(cols-1)*2];
                else
                    x_grid{end+1}=[idE,idE-(cols-1)*2,idE+2];
                end
            end

            idE=idE+1;
        end

    end
    i_grid=[i_grid,iList];
end

i_grid=i_grid(:);
end