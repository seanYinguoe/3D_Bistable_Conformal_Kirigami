function out = geobeam_discrete_2d(L,N,EI,EA,thetaA,thetaB,opts)
% Discrete geometrically-exact planar beam (Cosserat/Elastica)
% 离散几何精确二维梁（平面Cosserat/弹性线）— 变分能量最小化
%
% Inputs:
%   L      : length 长度
%   N      : number of segments 分段数 (nodes = N+1)
%   EI     : bending stiffness 弯曲刚度
%   EA     : axial stiffness 轴向刚度 (set 0 for inextensible 不可伸长≈0)
%   thetaA : end rotation at s=0 端部转角A (rad)
%   thetaB : end rotation at s=L 端部转角B (rad)
%   opts   : struct with fields (可选)
%       .Kt     tangent-consistency weight 切向一致性权重 (default 1e5*EI/L)
%       .Kbc    BC penalty 边界罚项 (default 1e6*EI/L)
%       .init   'straight'|'cubic' 初值 (default 'cubic')
%       .amp    plot amplification 放大倍数 (default 20)
%       .plot   true/false 是否绘图 (default true)
%       .solver 'auto'|'fminunc'|'gd' 求解器 (default 'auto')
%       .maxit  max iterations 最大迭代 (default 400)
%       .tol    gradient tol 梯度阈值 (default 1e-8)
%
% Output:
%   out struct: x,y,phi,energy,gradnorm,converged,history(可选)

if nargin<7, opts = struct; end
Kt   = getdef(opts,'Kt',  1e5*EI/max(L,1e-9));
Kbc  = getdef(opts,'Kbc', 1e6*EI/max(L,1e-9));
amp  = getdef(opts,'amp', 20);
mkp  = getdef(opts,'plot', true);
solver = getdef(opts,'solver','auto');
maxit = getdef(opts,'maxit',400);
tol   = getdef(opts,'tol',1e-8);
initType = getdef(opts,'init','cubic');

% ----- safety: allow degrees 自动识别角度制 -----
if max(abs([thetaA,thetaB])) > pi
    thetaA = deg2rad(thetaA); thetaB = deg2rad(thetaB);
    fprintf('[info] angles interpreted as degrees -> radians.\n');
end

h = L / N;

% ----- initial guess 初值 -----
switch lower(initType)
    case 'straight'
        x0 = linspace(0,L,N+1).'; y0 = zeros(N+1,1);
        phi0 = linspace(thetaA,thetaB,N+1).';
    otherwise % 'cubic' small-angle inspired
        x0 = linspace(0,L,N+1).'; 
        y0 = zeros(N+1,1);
        yy = ( (thetaA+thetaB)/L^2).*x0.^3 ...
           - ((2*thetaA+thetaB)/L).*x0.^2 + thetaA.*x0;
        y0 = yy;
        phi0 = atan2( gradient_s(yy,h), gradient_s(x0,h) + eps );
        phi0(1) = thetaA; phi0(end) = thetaB;
end

% pack variables: q = [x(1..N+1); y(1..N+1); phi(1..N+1)]
q0 = [x0; y0; phi0];

% choose solver 选择求解器
useFminunc = false;
if strcmpi(solver,'fminunc') || (strcmpi(solver,'auto') && exist('fminunc','file')==2)
    useFminunc = true;
end

if useFminunc
    Efun  = @(q) energy_and_grad(q,L,N,EI,EA,thetaA,thetaB,h,Kt,Kbc);
    optsOpt = optimoptions('fminunc','Algorithm','quasi-newton', ...
        'SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',maxit,...
        'OptimalityTolerance',tol,'StepTolerance',1e-12);
    [q, fval, exitflag, output] = fminunc(Efun, q0, optsOpt);
    converged = exitflag>0;
    hist = output;
else
    % simple gradient descent + backtracking 简易梯度下降（带回溯）
    q = q0; alpha0 = 1.0;
    [E, g] = energy_and_grad(q,L,N,EI,EA,thetaA,thetaB,h,Kt,Kbc);
    for it=1:maxit
        gnorm = norm(g);
        if gnorm < tol, break; end
        alpha = alpha0;
        c1 = 1e-4;
        % backtracking line-search 回溯直线搜索
        while true
            q_try = q - alpha*g;
            E_try = energy_only(q_try,L,N,EI,EA,thetaA,thetaB,h,Kt,Kbc);
            if E_try <= E - c1*alpha*(g.'*g), break; end
            alpha = alpha * 0.5;
            if alpha < 1e-12, break; end
        end
        q = q_try; E = E_try;
        [~, g] = energy_and_grad(q,L,N,EI,EA,thetaA,thetaB,h,Kt,Kbc);
        if it==1, alpha0 = alpha; else, alpha0 = 0.8*alpha + 0.2*alpha0; end
    end
    fval = E; converged = (norm(g)<tol);
    hist = struct('iterations',it,'final_grad',norm(g));
end

% unpack
[x,y,phi] = unpack_q(q,N);

% output
out = struct('x',x,'y',y,'phi',phi,'energy',fval, ...
             'converged',converged,'history',hist, ...
             'h',h,'N',N,'EI',EI,'EA',EA,'Kt',Kt,'Kbc',Kbc);

% plot
if mkp
    figure('Color','w'); hold on; grid on; box on;
    plot([0 L],[0 0],'k--','LineWidth',1);
    plot(x, amp*y,'b-','LineWidth',2);
    scatter([x(1) x(end)], amp*[y(1) y(end)], 40, 'r','filled');
    axis equal; xlim([0 L]); xlabel('x');
    ylabel(sprintf('y (amp=%g)',amp));
    title(sprintf('Discrete GE beam (N=%d), E=%.3e, conv=%d',N,fval,converged));
    legend({'baseline 基线','centerline 中心线','ends 端点'},'Location','best');
end

end

% ================= helpers =================

function [E, g] = energy_and_grad(q,L,N,EI,EA,thetaA,thetaB,h,Kt,Kbc)
% energy and analytic gradient 能量与解析梯度
[x,y,phi] = unpack_q(q,N);

% segment diffs
dx = x(2:end)-x(1:end-1);
dy = y(2:end)-y(1:end-1);
len = sqrt(dx.^2 + dy.^2) + 1e-15;

% stretch term 伸长项
lam = len / h;                             % lambda
Es = 0.5*EA*sum((lam-1).^2)*h;

% bending term 弯曲项 (wrap angle differences to [-pi,pi])
dphi = wrapToPi(phi(2:end)-phi(1:end-1));
kap  = dphi / h;
Eb   = 0.5*EI*sum(kap.^2)*h;

% tangent consistency 切向一致 (parallel transport penalty)
cx = dx - h*cos(phi(1:end-1));
cy = dy - h*sin(phi(1:end-1));
Et = 0.5*Kt*sum(cx.^2 + cy.^2);

% boundary penalties 边界罚项
Ebc = 0.5*Kbc*( (x(1)-0)^2 + (y(1)-0)^2 + (phi(1)-thetaA)^2 ...
              + (x(end)-L)^2 + (y(end)-0)^2 + (phi(end)-thetaB)^2 );

E = Es + Eb + Et + Ebc;

% ---- gradient 梯度（对 q = [x;y;phi]）----
% initialize
gx = zeros(N+1,1); gy=zeros(N+1,1); gphi = zeros(N+1,1);

% stretch grad
if EA~=0
    dE_dlen = EA*(lam-1)/h * h;  % = EA*(lam-1)
    % chain to x,y
    dlen_ddx = dx./len; dlen_ddy = dy./len;
    gx(1:end-1) = gx(1:end-1) - dE_dlen.*dlen_ddx;
    gx(2:end)   = gx(2:end)   + dE_dlen.*dlen_ddx;
    gy(1:end-1) = gy(1:end-1) - dE_dlen.*dlen_ddy;
    gy(2:end)   = gy(2:end)   + dE_dlen.*dlen_ddy;
end

% bending grad
dE_ddphi = EI*kap/h * h;  % = EI*kap/h * h = EI*kap/h? Careful:
% kap = dphi/h; Eb = 0.5*EI*sum((dphi/h).^2)*h = 0.5*EI/h * sum(dphi.^2)
% dEb/ddphi = (EI/h)*dphi
dE_ddphi = (EI/h)*dphi;
gphi(1:end-1) = gphi(1:end-1) - dE_ddphi;
gphi(2:end)   = gphi(2:end)   + dE_ddphi;

% tangent-consistency grad
% cx = dx - h*cos(phi_i), cy = dy - h*sin(phi_i)  with i=1..N
gx(1:end-1) = gx(1:end-1) - Kt*cx;
gx(2:end)   = gx(2:end)   + Kt*cx;
gy(1:end-1) = gy(1:end-1) - Kt*cy;
gy(2:end)   = gy(2:end)   + Kt*cy;
gphi(1:end-1) = gphi(1:end-1) - Kt*(-h*sin(phi(1:end-1)).*cx + h*cos(phi(1:end-1)).*cy);

% boundary gradients
gx(1)   = gx(1)   + Kbc*(x(1)-0);
gy(1)   = gy(1)   + Kbc*(y(1)-0);
gphi(1) = gphi(1) + Kbc*(phi(1)-thetaA);
gx(end)   = gx(end)   + Kbc*(x(end)-L);
gy(end)   = gy(end)   + Kbc*(y(end)-0);
gphi(end) = gphi(end) + Kbc*(phi(end)-thetaB);

g = [gx; gy; gphi];
end

function E = energy_only(q,L,N,EI,EA,thetaA,thetaB,h,Kt,Kbc)
E = energy_and_grad(q,L,N,EI,EA,thetaA,thetaB,h,Kt,Kbc);
E = E(1);
end

function [x,y,phi] = unpack_q(q,N)
m = N+1; x = q(1:m); y = q(m+1:2*m); phi = q(2*m+1:3*m);
end

function v = gradient_s(f,h)
% centered difference with end-caps 中心差分（端点前向/后向）
v = zeros(size(f));
v(2:end-1) = (f(3:end)-f(1:end-2))/(2*h);
v(1) = (f(2)-f(1))/h; v(end) = (f(end)-f(end-1))/h;
end

function v = getdef(S, name, defaultVal)
if isfield(S,name) && ~isempty(S.(name)), v = S.(name); else, v = defaultVal; end
end
