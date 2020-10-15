%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
function top88_heat(nelx,nely,volfrac,penal,rmin,ft)
tic
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-6;
%% PREPARE FINITE ELEMENT ANALYSIS
KE = [2/3 -1/6 -1/3 -1/6; -1/6 2/3 -1/6 -1/3; -1/3 -1/6 2/3 -1/6; -1/6 -1/3 -1/6 2/3];
ndof = (nely+1)*(nelx+1);
nele = nelx*nely;          % ALso use for nele
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nele,1);
edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
BC = 'Heat_Transfer';
stopping_criteria = 'change';
U = zeros(ndof,1);
emptyelts = [];
fullelts = [];
switch BC
    case 'Heat_Transfer'
        F = sparse(1:ndof,1,-0.01,ndof,1);
        fixeddofs = 1+(nely+1)*[floor((nelx+1)*0.4:1:(nelx+1)*0.6)]';
    otherwise
        error('BC string should be a valid entry: ''MBB'',''L-Shape'',''Short_Cantilever''')
end
mult = 1 - numel(emptyelts)/nele;       % Multiplier for volfrac in case of emptyelts
alldofs = [1:ndof]';
freedofs = setdiff(alldofs,fixeddofs);
%% Generate a *folder* and a prefix to save images *optimization history*:

rs=replace(num2str(rmin,'%3.2f'),'.','_');
folder_name=['Optimization_history_',BC,'_SIMP_','Volfrac',num2str(volfrac),'_nelx_',num2str(nelx),...
    '_nely_',num2str(nely),'_R_',rs];
image_prefix=[BC,'_nelx_',num2str(nelx),'nely_',num2str(nely),'_R_',rs];
mkdir(folder_name)
Path=[folder_name,'/'];
%% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
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
xmin = zeros(nele,1);
xmax = ones(nele,1);
xold1 = x(:);
xold2 = x(:);
low = ones(nele,1);
upp = ones(nele,1);
a0 = 1;
a = zeros(m,1);
c_MMA = 10000*ones(m,1);
d = zeros(m,1);
maxoutit = 2000;
kkttol  = 0.001;
changetol = 0.001;
kktnorm = kkttol+10;
change=1;
stop_cond = 1;
cvec = zeros(maxoutit,1);
vvec = cvec;
kvec = cvec;
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
print([Path,'density_000'],'-dpng')
%% Create a file to track iterations in case of HPC
f1 = fopen([image_prefix,'.txt'],'w+');
%% START ITERATION
while stop_cond
  iter = iter + 1;
  outeriter = outeriter + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+(E0-Emin)*xPhys(:)'.^penal),16*nele,1);
  K = sparse(iK,jK,sK); %K = (K+K')/2;
  U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx]);
  c = sum(sum((Emin+(E0-Emin)*xPhys.^penal).*ce));
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
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
  f0val = c;
  df0dx = dc(:);
  fval = sum(xPhys(:))/(volfrac*nele)-1;
  dfdx = dv(:)'/(volfrac*nele);
  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
      mmasub(m,nele,iter,xval,xmin,xmax,xold1,xold2, ...
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
        kktcheck(m,nele,xmma,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,c_MMA,d);
  % update the stopping criterion
    switch stopping_criteria
        case 'kktnorm'
            stop_cond=iter < maxoutit && kktnorm>kkttol;
            kvec(iter) = kktnorm;
        case 'change'
            stop_cond=iter < maxoutit &&change>changetol;
    end
    %% Greyness level
    gl = 4/nele*sum(xPhys(:).*(1-xPhys(:)));
  %% PRINT RESULTS
  fprintf(f1,' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',iter,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  if or(mod(iter,10) == 0,iter == 1)
      print([Path,'density_',num2str(iter,'%03d')],'-dpng')
      figure(2)
      subplot(2,1,1)
      plot(1:iter,cvec(1:iter),'bo','MarkerFaceColor','b')
      grid on
      hold on
      scatter(iter,c,'k','fill')
      hold off
      text(iter,c,['C =',num2str(c,'%4.2f'),' at iteration ', num2str(iter)],...
        'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
      xlabel('iter')
      ylabel('C')
      subplot(2,1,2)
      plot(1:iter,vvec(1:iter)*100,'ro','MarkerFaceColor','r')
      grid on
      hold on
      scatter(iter,mean(xPhys(:))*100,'k','fill')
      hold off
      text(iter,mean(xPhys(:))*100,['V = ',num2str(mean(xPhys(:))*100,'%4.2f'),'% at iteration ', num2str(iter)],...
      'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',24,'FontWeight','bold')
      xlabel('iter')
      ylabel('V [%]')
      print([Path,image_prefix,'_convergence'],'-dpng')
      % KKT PLOT
      figure(3)
      plot(1:iter,kvec(1:iter),'bo','MarkerFaceColor','b');
      grid on; hold on;
      scatter(iter,c,'k','fill')
      hold off
      xlabel('iter')
      ylabel('kktnorm')
      print([Path,image_prefix,'_kktnorm'],'-dpng')
  end
end
print([Path,image_prefix,'Density'],'-dpng');
fprintf(f1,'Elapsed Time: %10.3f',toc);
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

