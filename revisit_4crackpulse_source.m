clear all

dx = 0.025;
nx = 1601;
ny = 1601;
x = (0:nx-1)*dx-20;
y = (0:ny-1)*dx-20;

for i = 1:4
folder=['model',num2str(i)];
fid = fopen([folder,'/slip_x'],'r');
slip{i}=reshape(fread(fid,'float'),nx,ny);
fclose(fid);

fid = fopen([folder,'/tsoriginal'],'r');
ts0{i}=reshape(fread(fid,'float'),nx,ny);
fclose(fid);

fid = fopen([folder,'/tsfinal'],'r');
tse{i}=reshape(fread(fid,'float'),nx,ny);
fclose(fid);

strdrop{i}=ts0{i}-tse{i};
end


%% stress drop
G = 3.2038e+10;
subarea=ones(nx,ny)*dx*dx*1e6;

R=18e3;
A=15.2e3;
B=13.5e3;

xarray = repmat(x*1e3,ny,1);
yarray = repmat(y'*1e3,1,nx);

radius=sqrt(xarray.^2+yarray.^2);
radius(radius>=R)=R;

E12 = 24/7/pi*R*sqrt(1-radius.^2/R^2);

C={2.44,2.44,2.44,2.53};
for i =1:4
    M0 = G*slip{i}.*subarea;
    M0 = sum(M0(:));
    area=sum(subarea(slip{i}>1e-3),'all');
    
    fprintf('%e,%e\n',M0,area)
    
    tmp1=strdrop{i}.*slip{i}.*subarea;
    tmp2=slip{i}.*subarea;
    sigma1{i} = sum(tmp1(:))/sum(tmp2(:));
    if i<4
        sigma2{i}=M0*7/(16*R^3);
    else
        m=sqrt(1-B^2/A^2);
        [K,E] = ellipke(m);
        c1=4/(3*E+(E-B^2/A^2*K)/m^2);
        sigma2{i}=M0/(c1*area*B);
    end
    
    tmp1=strdrop{i}.*E12.*subarea;
    tmp2=E12.*subarea;
    sigma3{i} = sum(tmp1(:))/sum(tmp2(:));
    
    
    sigma4{i}=C{i}*M0/(area^(1.5));
    
end
fprintf('                      EC           GP           SP           AP\n')
fprintf('Method-1 (MPa): %12.5e,%12.5e,%12.5e,%12.5e\n',sigma1{1}/1e6,sigma1{2}/1e6,sigma1{3}/1e6,sigma1{4}/1e6)
fprintf('Method-2 (MPa): %12.5e,%12.5e,%12.5e,%12.5e\n',sigma2{1}/1e6,sigma2{2}/1e6,sigma2{3}/1e6,sigma2{4}/1e6)
fprintf('Method-3 (MPa): %12.5e,%12.5e,%12.5e\n',sigma3{1}/1e6,sigma3{2}/1e6,sigma3{3}/1e6)
fprintf('Method-4 (MPa): %12.5e,%12.5e,%12.5e,%12.5e\n',sigma4{1}/1e6,sigma4{2}/1e6,sigma4{3}/1e6,sigma4{4}/1e6)



%%





% figure(1)
% clf
% figure(2)
% clf
% 
% for i=1:4
% figure(1)
% subplot(2,2,i)
% colormap(parula)
% pcolor(x,y,slip{i}');
% shading flat;
% colorbar;
% axis equal
% xlim([-20,20])
% ylim([-20,20])
% 
% 
% 
% figure(2)
% subplot(2,2,i)
% colormap(turbo)
% pcolor(x,y,(-strdrop{i})'/1e6);
% shading flat;
% colorbar;
% caxis([-45,45])
% axis equal
% xlim([-20,20])
% ylim([-20,20])
% end



