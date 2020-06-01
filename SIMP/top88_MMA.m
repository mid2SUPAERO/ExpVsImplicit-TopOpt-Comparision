%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function top88_MMA(nelx,nely,volfrac,penal,rmin,ft,BC)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
ndof = 2*(nely+1)*(nelx+1);
n = nelx*nely;          % ALso use for nele
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
% BC = 'Short_Cantilever';
stopping_criteria = 'kktnorm';
U = sparse(ndof,1);
emptyelts = [];
fullelts = [];
mult = 1;       % Multiplier for volfrac in case of emptyelts
switch BC
    case 'MBB'
        F = sparse(2,1,-1,ndof,1);
        fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)])';
    case 'L_shape'
        fixedx = [1:2*(nely+1):ndof]';
        fixeddofs = union(fixedx,fixedx+1);
        ind = nely*[floor(nelx/2):1:nelx-1];
        add = [1:1:nely/2]';
        emptyelts = reshape(repmat(ind, size(add,1), 1) + repmat(add, 1, size(ind,2)),[],1);
        F = sparse( floor((emptyelts(end)-1)/nely+1)*(nely+1)*2 + mod(emptyelts(end),nely)*2+2, 1, -1, ndof,1);
        mult = 1 - numel(emptyelts)/n;
        Emin = 1e-6;
    case 'Short_Cantilever'
        F = sparse( 2*(nely+1)*(nelx)+floor((2*nely+1)/2)+2, 1, -1, ndof,1);
        fixeddofs = [1:1:2*(nely+1)]';
    case 'Compliant'
        in = 1; out = 2*(nely+1)*nelx+1;
        F = sparse([in;out],[1;2],[1; -1],ndof,2);
        U = sparse(ndof,2);
        fixedx = [2:2*(nely+1):ndof]';
        fixeddofs = union(fixedx,[2*(nely+1);2*(nely+1)-1]);
    otherwise
        error('BC string should be a valid entry: ''MBB'',''L_shape'',''Short_Cantilever''')
end

alldofs = [1:ndof]';
freedofs = setdiff(alldofs,fixeddofs);
%% Generate a *folder* and a prefix to save images *optimization history*:

rs=replace(num2str(rmin,'%3.2f'),'.','_');
vf = replace(num2str(volfrac,'%3.2f'),'.','_');
folder_name=['Optimization_history_',BC,'_SIMP_','Volfrac_',vf,'_nelx_',num2str(nelx),...
    '_nely_',num2str(nely),'_R_',rs,'_SC_',stopping_criteria];
image_prefix=[BC,'_SIMP_','Volfrac_',vf,'_nelx_',num2str(nelx),'nely_',num2str(nely),'_R_',rs];
mkdir(folder_name)
Path=[folder_name,'/'];
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
volfrac = mult*volfrac;
x = repmat(volfrac,nely,nelx);
xPhys = x;
iter = 0;
outeriter = 0;
%% INITIALIZE MMA OPTIMISER
m = 1;
xmin = zeros(n,1);
xmax = ones(n,1);
xold1 = x(:);
xold2 = x(:);
low = ones(n,1);
upp = ones(n,1);
a0 = 1;
a = zeros(m,1);
c_MMA = 10000*ones(m,1);
d = zeros(m,1);
maxoutit = 2000;
kkttol  = 0.0001;
changetol = 0.001;
kktnorm = kkttol+10;
change=1;
stop_cond = 1;
cvec = zeros(maxoutit,1);
vvec = cvec;
kvec = cvec;
figure(1); colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
print([Path,'density_000'],'-dpng');
jump = 5;
%% Write logs in file for tracking in case of HPC
f1 = fopen([image_prefix,'.txt'],'a+');
%% START ITERATION
while stop_cond
  iter = iter + 1;
  outeriter = outeriter + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin*xPhys(:)'+xPhys(:)'.^penal*(E0-Emin)),64*n,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  switch BC
      case 'Compliant'
          K(in,in) = K(in,in) + 0.1;
          K(out,out) = K(out,out) + 0.1;
  end
  U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  switch BC
      case 'Compliant'
          U1 = U(:,1); U2 = U(:,2);
          ce = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),[nely,nelx]);
          c  = U(out,1);
          dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce; %penal*(E0-Emin)*xPhys.^(penal-1).*ce;
      otherwise
          ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
          c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
          dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  end
  dv = ones(nely,nelx);
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  cvec(iter) = c;
  vvec(iter) = mean(xPhys(:));
  %% MMA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  xval = x(:);
  f0val = c*5;
  df0dx = dc(:)*5;
  fval = sum(xPhys(:))/(volfrac*n)-1;
  dfdx = dv(:)'/(volfrac*n);
  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
      mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
      f0val, df0dx, fval, dfdx, low, upp, a0, a, c_MMA,d);
  xmma(emptyelts) = Emin;
  xmma(fullelts) = E0;
  xnew = reshape(xmma,nely,nelx);
  xPhys(:) = (H*xnew(:))./Hs;
  xold2 = xold1(:);
  xold1 = x(:);
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  [residu,kktnorm,residumax] = ...
        kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,c_MMA,d);
  % update the stopping criterion
    switch stopping_criteria
        case 'kktnorm'
            stop_cond = (iter < maxoutit && kktnorm>kkttol);
            kvec(iter) = kktnorm;
        case 'change'
            stop_cond = (iter < maxoutit && change>changetol);
            kvec(iter) = change;
    end
  %% PRINT RESULTS
  fprintf(f1,' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f kktnrom: %7.4f\n',iter,full(c),mean(xPhys(:)),change,kktnorm);
  %% PLOT DENSITIES
  if iter < 101
      jump = 5;
  else
      jump = 20;
  end
  figure(1);
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  if or(mod(iter,jump) == 0,iter == 1)
      print([Path,'density_',num2str(iter,'%03d')],'-dpng');
      if iter > 20
          figure(2)
          subplot(2,1,1)
          plot(20:iter,cvec(20:iter),'bo','MarkerFaceColor','b')
          grid on
%           hold on
%           scatter(iter,c,'k','fill')
%           hold off
          text(iter,c,['C =',num2str(c,'%4.2f'),' at iteration ', num2str(iter)],...
            'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
          xlabel('iter')
          ylabel('C')
          subplot(2,1,2)
          plot(20:iter,vvec(20:iter)*100,'ro','MarkerFaceColor','r')
          grid on
%           hold on
%           scatter(iter,mean(xPhys(:))*100,'k','fill')
%           hold off
          text(iter,mean(xPhys(:))*100,['V = ',num2str(mean(xPhys(:))*100,'%4.2f'),'% at iteration ', num2str(iter)],...
          'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
          xlabel('iter')
          ylabel('V [%]')
          print([Path,image_prefix,'_convergence'],'-dpng')
          % KKT PLOT
          figure(3)
          plot(20:iter,kvec(20:iter),'bo','MarkerFaceColor','b');
          grid on;
%           hold on;
%           scatter(iter,c,'k','fill')
%           hold off
          xlabel('iter')
          switch stopping_criteria
          case 'kktnorm'
              ylabel('kktnorm')
              print([Path,image_prefix,'_kktnorm'],'-dpng')
          case 'change'
              ylabel('change')
              print([Path,image_prefix,'_change'],'-dpng')
          end
      end
  end
end
print('-f1',[Path,image_prefix,'Density'],'-dpng');
fclose(f1);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

