x = [0, 0.25, 0.5, 0.75, 1];
epsilon = [1, 1, 2, 2];
node_index = [1, 2, 3, 4;
              2, 3, 4, 5];
A = zeros(5,5);
element_A = zeros(2,2);
node = zeros(2,1);
dndx = zeros(2,1);

for element = 1:4
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

% A(1,:) = 0;
% A(5,:) = 0;
% A(5,5) = 1;

%B = zeros(5,1);
% format rat
% V = null(A, 'r');
% disp V