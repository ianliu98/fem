%% パラメータ
epsilon_0 = 8.854188 * 10^(-12); 	% 真空の誘電率
epsilon_r(1) = 1; 			% 材料1（真空）の比誘電率
epsilon_r(2) = 10; 		% 材料2（誘電体）の比誘電率
V_0 = 10.0; 			% 電極電位
%% メッシュ分割
xmin = 0;  xmax = 0.02;  noex = 20;	% [xmin xmax]をnoex分割する
ymin = 0;  ymax = 0.02;  noey = 20; 	% [ymin ymax]をnoey分割する
xaxis = linspace(xmin, xmax, noex + 1); 	% メッシュのx座標
yaxis = linspace(ymin, ymax, noey + 1); 	% メッシュのy座標
noe = 2 * noex * noey; 			% 三角形要素の数
nond = (noex + 1) * (noey + 1); 		% 節点の数
%% 材料の配置
material = ones(2, noex, noey); 	% 真空
material(:, 1:5, 1:2) = 2;		% 真空以外の材料の設定

%% 節点の座標・番号の設定
% 節点座標の設定
[x_node, y_node] = ndgrid( xaxis, yaxis );
% 各要素を構成する節点の番号（x, yの辞書式順）
element2node = zeros(3, 2, noex, noey);
for kx = 1:noex
    for ky = 1:noey
        % 左下の三角形
        element2node(1, 1, kx, ky) = kx + (ky - 1) * (noex + 1);
        element2node(2, 1, kx, ky) = kx + 1 + (ky - 1) * (noex + 1);
        element2node(3, 1, kx, ky) = kx + ky * (noex + 1);
        % 右上の三角形
        element2node(1, 2, kx, ky) = kx + 1 + (ky - 1) * (noex + 1);
        element2node(2, 2, kx, ky) = kx + 1 + ky * (noex + 1);
        element2node(3, 2, kx, ky) = kx + ky * (noex + 1);
    end
end
element2node = reshape(element2node, 3, noe);
 
%% 係数行列の計算
imat = ones(9, noe);  jmat = ones(9, noe);  mat = zeros(9, noe);
% 要素行列を計算する
for  i = 1:noe
    x = x_node( element2node(:, i) ); % 要素節点の x座標
    y = y_node( element2node(:, i) ); % 要素節点の y座標
    epsilon = epsilon_0 * epsilon_r( material(i) ); % 要素の誘電率
    a = [ x(2)*y(3)-x(3)*y(2);  x(3)*y(1)-x(1)*y(3);  x(1)*y(2)-x(2)*y(1) ];
    b = [ y(2)-y(3);  y(3)-y(1);  y(1)-y(2) ];
    c = [ x(3)-x(2);  x(1)-x(3);  x(2)-x(1) ];
    S = sum(a) / 2.0; % 要素の面積
    A_e = (b * b' + c * c') * epsilon / (4.0 * S) ; % 要素行列
    imat(:, i) = repmat( element2node(:, i), 3, 1 );
    jmat(:, i) = reshape( repmat( element2node(:, i), 1, 3 )' , 9, 1 );
    mat(:, i) = A_e(:);
end
% 係数行列は、要素行列の重ね合わせ
A = sparse(imat, jmat, mat);

%% 境界条件の設定
earth = 1:noex+1;  v_0 = [43:48 64:69]; % 電位 0, V_0 の条件を与える節点
dirichlet = [earth v_0];
A(dirichlet, :) = 0.0; % 固定境界の行を全て零に
A(dirichlet, dirichlet) = speye( length( dirichlet ) ); % 対角成分は 1 に
b = zeros(nond, 1);
b(v_0) = V_0;
 
%% 求解（内部ではUMFPACKが使用される）
x = A \ b;
 
%% 電位表示
[X, Y] = meshgrid( xaxis, yaxis );
Z = reshape(x, noex + 1, noey + 1)' ;
figure;  surf(X, Y, Z);  shading interp
figure;  contour(X, Y, Z) % 等高線

%% 電界表示
% 要素毎の電界
Ex = zeros(noe, 1); Ey = zeros(noe, 1);
X =  zeros(noe, 1); Y = zeros(noe, 1);
for i = 1:noe
    xe = x_node( element2node(:, i) ); % 要素節点の x座標
    ye = y_node( element2node(:, i) ); % 要素節点の y座標
    a = [ xe(2)*ye(3)-xe(3)*ye(2); xe(3)*ye(1)-xe(1)*ye(3); xe(1)*ye(2)-xe(2)*ye(1) ];
    b = [ ye(2)-ye(3); ye(3)-ye(1); ye(1)-ye(2) ];
    c = [ xe(3)-xe(2); xe(1)-xe(3); xe(2)-xe(1) ];
    S = sum(a) / 2.0; % 要素の面積
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
