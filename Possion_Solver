function [u_new, v_new, p_interior] = Possion(rho_cc, ustar, vstar, dx, dy, dt, ghostnum)

% Corrected Poisson projection for staggered grid.

[Ny, Nx] = size(rho_cc);
start = ghostnum + 1;
endy  = Ny - ghostnum;
endx  = Nx - ghostnum;
M = endy - start + 1;
N = endx - start + 1;
nCells = M * N;

% beta at cell centers
beta_cc = 1 ./ rho_cc;

% preallocate
I = zeros(5*nCells,1); J = I; S = I; ptr = 0;
b = zeros(nCells,1);

idx = @(jj,ii) (ii - start) * M + (jj - start) + 1;

% assemble A * p = b for interior cells
for j = start:endy
  for i = start:endx
    k = idx(j,i);

    % compute harmonic-average beta on faces
    small = 1e-14;
    beta_e = 2 * beta_cc(j,i)   * beta_cc(j,i+1) / (beta_cc(j,i)   + beta_cc(j,i+1)   + small);
    beta_w = 2 * beta_cc(j,i-1) * beta_cc(j,i)   / (beta_cc(j,i-1) + beta_cc(j,i)     + small);
    beta_n = 2 * beta_cc(j+1,i) * beta_cc(j,i)   / (beta_cc(j+1,i) + beta_cc(j,i)     + small);
    beta_s = 2 * beta_cc(j,i)   * beta_cc(j-1,i) / (beta_cc(j,i)   + beta_cc(j-1,i)   + small);

    aE = beta_e / dx^2;
    aW = beta_w / dx^2;
    aN = beta_n / dy^2;
    aS = beta_s / dy^2;

    diagA = (aE + aW + aN + aS);

    % RHS 
    divU = (ustar(j,i+1) - ustar(j,i)) / dx + (vstar(j+1,i) - vstar(j,i)) / dy;
    b(k) = -(1 / dt) * divU;

    % West neighbor
    if i-1 >= start
      ptr = ptr + 1; I(ptr)=k; J(ptr)=idx(j,i-1); S(ptr)=-aW;
    else
      % neighbor outside -> Neumann (zero normal) -> remove corresponding flux term
      % removing west flux is equivalent to NOT including aW* p_W term, and thus
      % diagonal should NOT include aW. Since diagA currently includes it, subtract it.
      diagA = diagA - aW;
    end

    % East neighbor
    if i+1 <= endx
      ptr = ptr + 1; I(ptr)=k; J(ptr)=idx(j,i+1); S(ptr)=-aE;
    else
      diagA = diagA - aE;
    end

    % South neighbor (j-1)
    if j-1 >= start
      ptr = ptr + 1; I(ptr)=k; J(ptr)=idx(j-1,i); S(ptr)=-aS;
    else
      diagA = diagA - aS;
    end

    % North neighbor (j+1)
    if j+1 <= endy
      ptr = ptr + 1; I(ptr)=k; J(ptr)=idx(j+1,i); S(ptr)=-aN;
    else
      diagA = diagA - aN;
    end

    ptr = ptr + 1; I(ptr)=k; J(ptr)=k; S(ptr)=diagA;
  end
end

% trim
I = I(1:ptr); J = J(1:ptr); S = S(1:ptr);
A = sparse(I,J,S,nCells,nCells);

% enforce reference p = 0 at one interior cell
% choose the central interior cell as reference
ref_i = round((start+endx)/2);
ref_j = round((start+endy)/2);
ref_k = idx(ref_j, ref_i);
% replace ref_k row by identity
A(ref_k, :) = 0;
A(ref_k, ref_k) = 1;
b(ref_k) = 0;

% Solve
pvec = A \ b;

% reshape to full pcc with ghosts
pcc = zeros(Ny, Nx);
p_interior = reshape(pvec, M, N);
pcc(start:endy, start:endx) = p_interior;


% compute face-centered rho for correction 
rho_u = ones(Ny+1, Nx+1) * eps;
rho_v = ones(Ny+1, Nx+1) * eps;
for j = start:endy+1
  for i = start:endx+1
    % use neighboring cell-centered rho
    rhoL = rho_cc(j, i-1); rhoR = rho_cc(j, i);
    % harmonic-like: rho_face = 2/(1/rhoL + 1/rhoR) = 2 * rhoL * rhoR / (rhoL + rhoR)
    rho_u(j,i) = 2 * rhoL * rhoR / (rhoL + rhoR + 1e-14);
    rhoB = rho_cc(j-1, i); rhoT = rho_cc(j, i);
    rho_v(j,i) = 2 * rhoB * rhoT / (rhoB + rhoT + 1e-14);
  end
end

% Correct velocities
u_new = ustar;
v_new = vstar;
for j = start:endy+1
  for i = start:endx+1
    dpdx = ( pcc(j, i) - pcc(j, i-1) ) / dx;
    u_new(j,i) = ustar(j,i) - dt * dpdx / rho_u(j,i);
    dpdy = ( pcc(j, i) - pcc(j-1, i) ) / dy;
    v_new(j,i) = vstar(j,i) - dt * dpdy / rho_v(j,i);
  end
end

% Boundary enforcement 
% set wall face velocities to zero (no-slip)

%left
for j=1:Ny+1
    for i=1:ghostnum
        u_new(j,i)=-u_new(j,2*ghostnum+2-i);
        v_new(j,i)=-v_new(j,2*ghostnum+1-i);
    end
end

%right
for j=1:Ny+1
    for i=endx+1:Nx+1
        u_new(j,i)=-u_new(j,2*endx+2-i);
    end
end

for j=1:Ny+1
    for i=endx+1:Nx+1
        v_new(j,i)=-v_new(j,2*endx+1-i);
    end
end

%up
for j=1:ghostnum
    for i=1:Nx+1
        u_new(j,i)=-u_new(2*ghostnum+1-j,i);
        v_new(j,i)=-v_new(2*ghostnum+2-j,i);
    end
end

%down
for j=endy+1:Ny+1
    for i=1:Nx+1
        u_new(j,i)=-u_new(2*endy+1-j,i);
        
    end
end

for j=endy+1:Ny+1
    for i=1:Nx+1
        v_new(j,i)=-v_new(2*endy+2-j,i);
    end
end


u_new(:, start)   = 0;
u_new(:, endx+1)  = 0;
v_new(start, :)   = 0;
v_new(endy+1, :)  = 0;

% return interior pressure (for diagnostic)
p_interior = pcc(start:endy, start:endx);

end
