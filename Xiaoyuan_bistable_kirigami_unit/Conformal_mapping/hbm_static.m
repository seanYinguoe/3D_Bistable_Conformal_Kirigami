function [Etotal, XY, out] = hbm_static(A0,B0,A1,B1,thetaA,thetaB,N,Kb,Ks,varargin)
% HBM_STATIC - Hencky bar-chain static large-deformation solver (静态大变形)
% Discretize a 2D Euler-Bernoulli-like beam as N rigid bars + (N-1) hinges.
% Energy = bending (hinges) + axial (bars) + penalties for end position/orientation.
%
% Inputs:
%   A0,B0   : [1x2] undeformed end points (未变形端点)
%   A1,B1   : [1x2] deformed end points   (变形端点) - target/goal
%   thetaA  : scalar, prescribed rotation at A (A端转角, rad)
%   thetaB  : scalar, prescribed rotation at B (B端转角, rad)
%   N       : integer, number of bars (杆段数)
%   Kb      : bending stiffness per hinge (每个铰的弯曲刚度)
%   Ks      : axial stiffness per bar     (每段的轴向刚度)
%
% Name-Value options:
%   'w_pos'      : penalty for end-position mismatch (端点位置罚权重), default 1e6
%   'w_ang'      : penalty for end-angle mismatch    (端点转角罚权重), default 1e6
%   'use_fminunc': use fminunc if available (若有工具箱优先用), default true
%   'x0'         : initial guess struct with fields psi (N-1), eps (N) (初值)
%
% Outputs:
%   Etotal : total energy at optimum (总能量)
%   XY     : (N+1)x2 node coordinates (节点坐标)
%   out    : struct with fields:
%            .Ebend, .Eaxial, .Epen_pos, .Epen_ang, .psi, .eps, .beta
%            .info (optimizer info), .opts (used options)

opts = struct('w_pos',1e6, 'w_ang',1e6, 'use_fminunc',true, 'x0',[]);
opts = parseOpts(opts, varargin{:});

% Geometry / 几何
A0 = A0(:).'; B0 = B0(:).';
A1 = A1(:).'; B1 = B1(:).';
L0 = norm(B0 - A0);         % undeformed length (未变形长度)
a0 = L0 / N;                % bar rest length (每段原长)
dTarget = (B1 - A1);        % target end-to-end vector in deformed config (目标端到端向量)

% Unknowns:
%   psi(1..N-1) : hinge angles (相邻杆之间的转角, rad)
%   eps(1..N)   : axial strain per bar (每段轴向应变), so current length = a0*(1+eps_i)
nPsi  = N-1;
nEps  = N;

% Initial guess / 初值
if isempty(opts.x0)
    psi0 = zeros(nPsi,1);                  % start straight (初始直杆)
    eps0 = (norm(dTarget)/L0 - 1)*ones(nEps,1); % uniform stretch to meet length scale
else
    % allow partial fields
    psi0 = getfield_or(opts.x0,'psi', zeros(nPsi,1));
    eps0 = getfield_or(opts.x0,'eps', zeros(nEps,1));
    psi0 = psi0(:); eps0 = eps0(:);
end
x0 = [psi0; eps0];

% Objective / 目标函数
obj = @(x) total_energy(x, N, a0, Kb, Ks, thetaA, thetaB, A1, dTarget, opts.w_pos, opts.w_ang);

% Solve / 求解
use_fminunc = opts.use_fminunc && exist('fminunc','file')==2;
if use_fminunc
    try
        o = optimoptions('fminunc','Algorithm','quasi-newton','Display','off',...
                         'MaxFunEvals',2e5,'MaxIter',5e3);
        [x, fval, exitflag, output] = fminunc(obj, x0, o);
        info.solver = 'fminunc'; info.exitflag = exitflag; info.output = output;
    catch
        [x, fval, exitflag, output] = fminsearch(obj, x0, optimset('Display','off','MaxFunEvals',2e5,'MaxIter',5e3));
        info.solver = 'fminsearch(fallback)'; info.exitflag = exitflag; info.output = output;
    end
else
    [x, fval, exitflag, output] = fminsearch(obj, x0, optimset('Display','off','MaxFunEvals',2e5,'MaxIter',5e3));
    info.solver = 'fminsearch'; info.exitflag = exitflag; info.output = output;
end

% Rebuild state / 重建结果
psi = x(1:nPsi);
eps = x(nPsi+1:end);
[~, XY, parts, beta] = total_energy(x, N, a0, Kb, Ks, thetaA, thetaB, A1, dTarget, opts.w_pos, opts.w_ang);

Etotal   = fval;
out.Ebend    = parts.Ebend;
out.Eaxial   = parts.Eaxial;
out.Epen_pos = parts.Epen_pos;
out.Epen_ang = parts.Epen_ang;
out.psi  = psi(:);
out.eps  = eps(:);
out.beta = beta(:);
out.info = info;
out.opts = opts;

end

% ------------ helpers / 辅助函数 ----------------

function [E, XY, parts, beta] = total_energy(x, N, a0, Kb, Ks, thetaA, thetaB, A1, dTarget, w_pos, w_ang)
% Unpack / 解包
nPsi = N-1;
psi  = x(1:nPsi);
eps  = x(nPsi+1:end);

% Absolute bar orientations from left end A (绝对杆方向角)
% beta(1) = thetaA; beta(i) = thetaA + sum_{k=1}^{i-1} psi(k)
beta = zeros(N,1);
beta(1) = thetaA;
if nPsi > 0
    beta(2:end) = thetaA + cumsum(psi(:));
end

% Bar lengths (当前每段长度)
li = a0*(1 + eps(:));

% Node coordinates (从 A1 开始累加, 得到变形后的节点)
XY = zeros(N+1,2);
XY(1,:) = A1;
for i = 1:N
    di = [cos(beta(i)), sin(beta(i))];
    XY(i+1,:) = XY(i,:) + li(i).*di;
end

% End-to-end closure error (端到端误差)
dNow = XY(end,:) - XY(1,:);
rPos = dNow - dTarget;                 % position mismatch (位置不匹配)
Epen_pos = 0.5*w_pos*dot(rPos, rPos);  % position penalty energy

% End orientation mismatch at B (末端转角不匹配)
% Target absolute orientation at the last bar is thetaB
dang = wrapToPi(beta(end) - thetaB);   % keep in (-pi,pi]
Epen_ang = 0.5*w_ang*(dang^2);         % angle penalty energy

% Bending energy (弯曲能, 每个内部铰)
Ebend = 0.5*Kb * sum(psi.^2);

% Axial energy (轴向能, 每段)
Eaxial = 0.5*Ks * sum( (li - a0).^2 ); % = 0.5*(Ks)*(a0^2)*eps.^2

% Total energy
E = Ebend + Eaxial + Epen_pos + Epen_ang;

% Pack parts
parts = struct('Ebend',Ebend,'Eaxial',Eaxial,'Epen_pos',Epen_pos,'Epen_ang',Epen_ang);
end

function val = getfield_or(S, name, defaultVal)
if isfield(S, name) && ~isempty(S.(name))
    val = S.(name);
else
    val = defaultVal;
end
end

function S = parseOpts(S, varargin)
if mod(numel(varargin),2)~=0
    error('Name-Value pairs expected.');
end
for k = 1:2:numel(varargin)
    nm = varargin{k};
    vl = varargin{k+1};
    if isfield(S, nm)
        S.(nm) = vl;
    else
        S.(nm) = vl; % allow extra
    end
end
end

function a = wrapToPi(a)
% Map angles to (-pi, pi] (角度映射)
a = mod(a + pi, 2*pi) - pi;
end
