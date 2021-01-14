epsilon = 10; % permittivity
V_0 = 1;
node_index = [1, 2, 3, 1, 1, 4, 6, 1, 1,  8,  9,  1;
              2, 3, 4, 4, 4, 6, 7, 7, 8,  9,  10, 10;
              5, 5, 5, 5, 7, 7, 8, 8, 10, 10, 11, 11]; % node
% x position         
x_1  = 0.5/sqrt(2);
x_2 = 0.5;
x_3 = 1.0;
x_4 = 1.0;
x_5 = (1 - 0.5*x_1)/(1.5 - x_1);
x_6 = 1.0;
x_7 = 0.5 * (x_1 + 1);
x_8 = x_1;
x_9 = 0;
x_10 = x_1/0.5 * ((1 - 0.5*x_1)/(1.5 - x_1) - 0.5);
x_11 = 0;
% y position
y_1 = 0.5/sqrt(2);
y_2 = 0;
y_3 = 0;
y_4 = y_1;
y_5 = y_1/0.5 * ((1 - 0.5*y_1)/(1.5 - y_1) - 0.5);
y_6 = 1.0;
y_7 = 0.5 * (y_1 + 1);
y_8 = 1.0;
y_9 = 1.0;
y_10 = (1 - 0.5*y_1)/(1.5 - y_1);
y_11 = 0.5;

x_node = [x_1, x_2, x_3, x_1, x_1, x_4, x_6, x_1, x_1,  x_8,  x_9,  x_1;
          x_2, x_3, x_4, x_4, x_4, x_6, x_7, x_7, x_8,  x_9,  x_10, x_10;
          x_5, x_5, x_5, x_5, x_7, x_7, x_8, x_8, x_10, x_10, x_11, x_11];
y_node = [y_1, y_2, y_3, y_1, y_1, y_4, y_6, y_1, y_1,  y_8,  y_9,  y_1;
          y_2, y_3, y_4, y_4, y_4, y_6, y_7, y_7, y_8,  y_9,  y_10, y_10;
          y_5, y_5, y_5, y_5, y_7, y_7, y_8, y_8, y_10, y_10, y_11, y_11];
    
noe = 12; % element number
nond = 11; % node number
imat = ones(9, noe);  jmat = ones(9, noe);  mat = zeros(9, noe); % store
% calculate coefficient matrix
for  i = 1:noe
    x = x_node(:,i); 
    y = y_node(:,i); 
    a = [ x(2)*y(3)-x(3)*y(2);  x(3)*y(1)-x(1)*y(3);  x(1)*y(2)-x(2)*y(1) ];
    b = [ y(2)-y(3);  y(3)-y(1);  y(1)-y(2) ];
    c = [ x(3)-x(2);  x(1)-x(3);  x(2)-x(1) ];
    S = sum(a) / 2.0; 
    A_e = (b * b' + c * c') * epsilon / (4.0 * S) ;
    imat(:, i) = repmat( node_index(:, i), 3, 1 );
    jmat(:, i) = reshape( repmat( node_index(:, i), 1, 3 )' , 9, 1 );
    mat(:, i) = A_e(:);
end
A = sparse(imat, jmat, mat); % sparse matrix
% boundary condition
earth = [1, 2, 11];
v_0 = [3, 4, 6, 8, 9];
dirichlet = [earth v_0];
A(dirichlet, :) = 0.0;
A(dirichlet, dirichlet) = speye( length( dirichlet ) );
b = zeros(nond, 1);
b(v_0) = V_0;
result = A \ b;

n1 = [x_1,y_1,rr(1)];
n2 = [x_2,y_2,rr(2)];
n3 = [x_3,y_3,rr(3)];
n4 = [x_4,y_4,rr(4)];
n5 = [x_5,y_5,rr(5)];
n6 = [x_6,y_6,rr(6)];
n7 = [x_7,y_7,rr(7)];
n8 = [x_8,y_8,rr(8)];
n9 = [x_9,y_9,rr(9)];
n10 = [x_10,y_10,rr(10)];
n11 = [x_11,y_11,rr(11)];

figure,
%axis([-0.1 1.1 -0.1 1.1 -1 1])
grid on
scatter3(x_1,y_1,rr(1),'*b')
hold on
scatter3(x_2,y_2,rr(2),'*b')
scatter3(x_3,y_3,rr(3),'*b')
scatter3(x_4,y_4,rr(4),'*b')
scatter3(x_5,y_5,rr(5),'*b')
scatter3(x_6,y_6,rr(6),'*b')
scatter3(x_7,y_7,rr(7),'*b')
scatter3(x_8,y_8,rr(8),'*b')
scatter3(x_9,y_9,rr(9),'*b')
scatter3(x_10,y_10,rr(10),'*b')
scatter3(x_11,y_11,rr(11),'*b')
plot3([n1(1), n2(1)], [n1(2), n2(2)], [n1(3), n2(3)], 'b')
plot3([n1(1), n11(1)], [n1(2), n11(2)], [n1(3), n11(3)], 'b')
plot3([n1(1), n10(1)], [n1(2), n10(2)], [n1(3), n10(3)], 'b')
plot3([n11(1), n10(1)], [n11(2), n10(2)], [n11(3), n10(3)], 'b')
plot3([n1(1), n7(1)], [n1(2), n7(2)], [n1(3), n7(3)], 'b')
plot3([n1(1), n5(1)], [n1(2), n5(2)], [n1(3), n5(3)], 'b')
plot3([n5(1), n2(1)], [n5(2), n2(2)], [n5(3), n2(3)], 'b')
plot3([n8(1), n7(1)], [n8(2), n7(2)], [n8(3), n7(3)], 'b')
plot3([n6(1), n7(1)], [n6(2), n7(2)], [n6(3), n7(3)], 'b')
plot3([n4(1), n7(1)], [n4(2), n7(2)], [n4(3), n7(3)], 'b')
plot3([n5(1), n4(1)], [n5(2), n4(2)], [n5(3), n4(3)], 'b')
plot3([n5(1), n3(1)], [n5(2), n3(2)], [n5(3), n3(3)], 'b')
plot3([n8(1), n10(1)], [n8(2), n10(2)], [n8(3), n10(3)], 'b')
plot3([n9(1), n10(1)], [n9(2), n10(2)], [n9(3), n10(3)], 'b')
plot3([n3(1), n4(1)], [n3(2), n4(2)], [n3(3), n4(3)], 'b')
plot3([n6(1), n4(1)], [n6(2), n4(2)], [n6(3), n4(3)], 'b')
plot3([n6(1), n8(1)], [n6(2), n8(2)], [n6(3), n8(3)], 'b')
plot3([n9(1), n8(1)], [n9(2), n8(2)], [n9(3), n8(3)], 'b')

title("electric potential")
xlabel('x')
ylabel('y')
zlabel('V(x,y)')
