%% �p�����[�^
epsilon_0 = 8.854188 * 10^(-12); 	% �^��̗U�d��
epsilon_r(1) = 1; 			% �ޗ�1�i�^��j�̔�U�d��
epsilon_r(2) = 10; 		% �ޗ�2�i�U�d�́j�̔�U�d��
V_0 = 10.0; 			% �d�ɓd��
%% ���b�V������
xmin = 0;  xmax = 0.02;  noex = 20;	% [xmin xmax]��noex��������
ymin = 0;  ymax = 0.02;  noey = 20; 	% [ymin ymax]��noey��������
xaxis = linspace(xmin, xmax, noex + 1); 	% ���b�V����x���W
yaxis = linspace(ymin, ymax, noey + 1); 	% ���b�V����y���W
noe = 2 * noex * noey; 			% �O�p�`�v�f�̐�
nond = (noex + 1) * (noey + 1); 		% �ߓ_�̐�
%% �ޗ��̔z�u
material = ones(2, noex, noey); 	% �^��
material(:, 1:5, 1:2) = 2;		% �^��ȊO�̍ޗ��̐ݒ�

%% �ߓ_�̍��W�E�ԍ��̐ݒ�
% �ߓ_���W�̐ݒ�
[x_node, y_node] = ndgrid( xaxis, yaxis );
% �e�v�f���\������ߓ_�̔ԍ��ix, y�̎��������j
element2node = zeros(3, 2, noex, noey);
for kx = 1:noex
    for ky = 1:noey
        % �����̎O�p�`
        element2node(1, 1, kx, ky) = kx + (ky - 1) * (noex + 1);
        element2node(2, 1, kx, ky) = kx + 1 + (ky - 1) * (noex + 1);
        element2node(3, 1, kx, ky) = kx + ky * (noex + 1);
        % �E��̎O�p�`
        element2node(1, 2, kx, ky) = kx + 1 + (ky - 1) * (noex + 1);
        element2node(2, 2, kx, ky) = kx + 1 + ky * (noex + 1);
        element2node(3, 2, kx, ky) = kx + ky * (noex + 1);
    end
end
element2node = reshape(element2node, 3, noe);
 
%% �W���s��̌v�Z
imat = ones(9, noe);  jmat = ones(9, noe);  mat = zeros(9, noe);
% �v�f�s����v�Z����
for  i = 1:noe
    x = x_node( element2node(:, i) ); % �v�f�ߓ_�� x���W
    y = y_node( element2node(:, i) ); % �v�f�ߓ_�� y���W
    epsilon = epsilon_0 * epsilon_r( material(i) ); % �v�f�̗U�d��
    a = [ x(2)*y(3)-x(3)*y(2);  x(3)*y(1)-x(1)*y(3);  x(1)*y(2)-x(2)*y(1) ];
    b = [ y(2)-y(3);  y(3)-y(1);  y(1)-y(2) ];
    c = [ x(3)-x(2);  x(1)-x(3);  x(2)-x(1) ];
    S = sum(a) / 2.0; % �v�f�̖ʐ�
    A_e = (b * b' + c * c') * epsilon / (4.0 * S) ; % �v�f�s��
    imat(:, i) = repmat( element2node(:, i), 3, 1 );
    jmat(:, i) = reshape( repmat( element2node(:, i), 1, 3 )' , 9, 1 );
    mat(:, i) = A_e(:);
end
% �W���s��́A�v�f�s��̏d�ˍ��킹
A = sparse(imat, jmat, mat);

%% ���E�����̐ݒ�
earth = 1:noex+1;  v_0 = [43:48 64:69]; % �d�� 0, V_0 �̏�����^����ߓ_
dirichlet = [earth v_0];
A(dirichlet, :) = 0.0; % �Œ苫�E�̍s��S�ė��
A(dirichlet, dirichlet) = speye( length( dirichlet ) ); % �Ίp������ 1 ��
b = zeros(nond, 1);
b(v_0) = V_0;
 
%% �����i�����ł�UMFPACK���g�p�����j
x = A \ b;
 
%% �d�ʕ\��
[X, Y] = meshgrid( xaxis, yaxis );
Z = reshape(x, noex + 1, noey + 1)' ;
figure;  surf(X, Y, Z);  shading interp
figure;  contour(X, Y, Z) % ������

%% �d�E�\��
% �v�f���̓d�E
Ex = zeros(noe, 1); Ey = zeros(noe, 1);
X =  zeros(noe, 1); Y = zeros(noe, 1);
for i = 1:noe
    xe = x_node( element2node(:, i) ); % �v�f�ߓ_�� x���W
    ye = y_node( element2node(:, i) ); % �v�f�ߓ_�� y���W
    a = [ xe(2)*ye(3)-xe(3)*ye(2); xe(3)*ye(1)-xe(1)*ye(3); xe(1)*ye(2)-xe(2)*ye(1) ];
    b = [ ye(2)-ye(3); ye(3)-ye(1); ye(1)-ye(2) ];
    c = [ xe(3)-xe(2); xe(1)-xe(3); xe(2)-xe(1) ];
    S = sum(a) / 2.0; % �v�f�̖ʐ�
    X(i) = sum(xe)/3.0;
    Y(i) = sum(ye)/3.0;
    Ex(i) = - x( element2node(:,i) )' * b / (2*S);
    Ey(i) = - x( element2node(:,i) )' * c / (2*S);
end
domain=[1:16 41:56 81:96 121:136 141:156];
figure; quiver(X(domain), Y(domain), Ex(domain), Ey(domain), 0.5)
set(gca, 'FontSize', 16)
hold on
plot([0 0.005 0.005 0],[0.002 0.002 0.003 0.003])
axis equal
axis([ 0.000, 0.008, - 0.001, 0.0045 ])
