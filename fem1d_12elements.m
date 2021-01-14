x = [0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1]; % partition
epsilon = [1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1]; % permittivity
node_index = [1, 2, 3, 4, 5, 6, 7, 8, 9,  10, 11, 12; % node
              2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
A = zeros(13,13);
element_A = zeros(2,2);
node = zeros(2,1);
dndx = zeros(2,1);

% calculate A matrix
for element = 1:12
   node(1) = node_index(1, element);
   node(2) = node_index(2, element);
   width = x(node(2)) - x(node(1));
   dndx(1) = -1/width;
   dndx(2) = 1/width;
   
   for i = 1:2
       for j= 1:2
           element_A(i,j) = epsilon(element) * dndx(i) * dndx(j) * width;
       end
   end
   
   for i = 1:2
       for j = 1:2
           A(node(i),node(j)) = A(node(i), node(j)) + element_A(i,j);
       end
   end
end

% calculate result
A_seg = A(2:12,2:12);
B = zeros(11,1);
B(11,1) = -A(12,13);
V = zeros(13,1);
V(13) = 1;
V(2:12) = A_seg\B;

% plot
figure,
plot(V, x)
title("Electric potential V")
axis([ 0.000, 1.100, 0.000, 1.100 ])
grid on
ylabel("V(x)")
xlabel("x")