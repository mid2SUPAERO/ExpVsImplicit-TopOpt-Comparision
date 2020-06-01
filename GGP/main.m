y= reshape(repmat([100,50,50],12,1),[],1);
volfrac = repmat(reshape(repmat([0.8 0.5 0.4 0.2],3,1),[],1),3,1);
shp = reshape(repmat([1 2 3],12,1),[],1);
mtd = reshape(repmat([1;2;3],1,12),[],1);
shpn = {'L-shape';'Short_Cantilever';'Compliant'};
mtdn = {'MMC';'MNA';'GP'};
fpar = [y volfrac shp mtd];
%%
parpool('local',12);
parfor i=[1:12]
    temp = fpar(i,:);
    GGP_main(100,temp(1),temp(2),shpn{temp(3)},mtdn{temp(4)})
end
% parfor i=[1:12]
%     temp = fpar(i,:);
%     GGP_heat(100,temp(1),temp(2),mtdn{temp(4)})
% end