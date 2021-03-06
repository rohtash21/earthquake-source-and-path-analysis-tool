% -------------------------------------------------------------------------
%           PROGRAM TO ESTIMATE fmax, p and kappa
%                   
% -------------------------------------------------------------------------
%       LOADING THE INPUT WAVEFORM DATA IN K-net FORMAT
%--------------------------------------------------------------------------
clc;clear all
for fff=1:3   
fileToRead1 = [num2str(fff) '.txt' ];
textfilename=fileToRead1;
fid = fopen(textfilename,'r');
data=textscan(fid,'%f','delimiter','\b','headerlines', 17);
data{1};
DATA=[data{1}*4.7648e+05];%4.7648e-04
fclose(fid);
%
z=DATA(:,1);
z=z-mean(z);
VX=z;
% -------------------------------------------------------------------------
SR=100;%input('Enter the sampling rate, SR: ');
SP=1/SR; % Sampling Period
N=1024;%input('Enter the number of samples,N : ');
% -------------------------------------------------------------------------
%--------------------------------------------------------------------------
scrsz = get(0,'ScreenSize');%
figure ('position',[scrsz]);
plot(z);hold on
CT=ginput(2);hold off;close
Mpc=round(CT(1));
Msc=round(CT(2));
SS=round(Mpc-1.5*(Msc-Mpc));
if SS<1
    SS=1;
else
        SS=SS; 
end
SE=round((Msc+2.4*(Msc-Mpc))+2*N);
Lz=length(z);
if SE>Lz
    SE=Lz;
else
        SE=SE;
end
VX=z(SS:SE);
%--------------------------------------------------------------------------
scrsz = get(0,'ScreenSize');% To assign the size of plotting window
figure ('position',[scrsz]);
subplot(2,1,1);plot(VX);title('Select P and S Onsets');grid on;
%-------------------------------------------------------------------------
MM=ginput(2);hold off;close
Mp=MM(1);
Ms=MM(2);
M=Ms;
tsp=(Ms-Mp)/SR;%30.12;
tsr=tsp+1.38*tsp;
gr=tsp*8*100000;
ln=length(VX);
t=(0:SP:(ln-1)*SP)';
tmin=0;tmax=max(t);
%--------------------------------------------------------------------------
%                       WINDOW FUNCTION 
%--------------------------------------------------------------------------
w=kaiser(N);%     COSINE TAPER  (20%)tukeywin(N,0.2);%
% for j=1:3
wdata(:,1)=VX(M:M+N-1,1).*w;
% end
% -------------------------------------------------------------------------
MX=max(abs(VX));
scrsz = get(0,'ScreenSize');% To assign the size of plotting window
figure ('position',[scrsz]);
subplot(3,1,1);plot(t,VX);title(['Accleration Time History:  ',num2str(fff)]);
xlabel('time(s)');ylabel('cm/s^2');axis([0 tmax -MX MX])
grid on;hold on

lt=max(MX);ht=2*lt;
rectangle('Position', [M*SP -lt N*SP ht],'FaceColor','g'...
               ,'EdgeColor','red','LineWidth',1);
plot(t(M+1:M+N),VX(M+1:M+N),'r');
% -------------------------------------------------------------------------
saveas(gcf,['PlotTH' num2str(fff)],'fig');
saveas(gcf,['plotTH' num2str(fff)],'bmp');hold off;close
% -------------------------------------------------------------------------
n=N/2;  
f=(1:n)*(SR/N);
tsp=(M-MM)/100;%30.12;
tsr=tsp+1.38*tsp;
gr=tsp*8*100000;
datav(:,1)=f;
%for i=1:3
    df=wdata(:,1);
    b=fft(df);
    b=abs(b(1:n))./n;
% %     q=5.4.*(1+f/0.3).^3.5./(f/0.3).^2;%82.*f.^1.12;%
% %     atnn=exp(-pi*f*tsr./q);%a=exp(nu./de)
    bb=b';%.*atnn;
    datav(:,2)=bb./(2*pi*f);
%   end
fln = [ 'SPECVEL',num2str(fff) '.txt' ];
dlmwrite(fln,[ f bb ],'delimiter','\t','precision',6)
%--------------------------------------------------------------------------
%        PROGRAM TO ESTIMATE SOURCE PARAMETERS
%--------------------------------------------------------------------------
fv=datav(:,1);
vel=datav(:,2);
displ=vel./(2*pi*fv).^1;acc=vel.*(2*pi*fv).^1;
loglog(fv,acc,'m');grid on;title(['Accleration Spectrum:  ',num2str(fff)]);
xlabel('frequency (Hz)'); ylabel( 'Acceleration Spectrum');grid on;hold on
% =========================================================================
%                ESTMATION OF Omega0, fc and fmax
% -------------------------------------------------------------------------
spec=ginput(4);
fc=spec(2);
fmax=spec(3);
omg1=spec(2,2);
omg2=spec(3,2);
avacc=(omg1+omg2)/2;
omg=avacc/(2*pi*fc).^2;
omg=omg+0.5*omg;
% -------------------------------------------------------------------------
%                   Estimation of slope above fmax
%--------------------------------------------------------------------------
px1=spec(3,1);
py1=spec(3,2);
px2=spec(4,1);
py2=spec(4,2);
slop=log10(py1/py2)/(log10(px1/px2));
p=abs(slop);
%--------------------------------------------------------------------------
%       Evaluation of kappa according to least error fit
%--------------------------------------------------------------------------
    kappa=0.001:0.001:0.08;
    lkappa=length(kappa);
    for i=1:lkappa
        kappa1=kappa(i);
    outkappa(1,i)=kappa1;
    spektraf=omg./((1.+(fv./fc).^2).*exp(pi*fv*kappa1));
    outspek(:,i)=spektraf;
    defr(i)=mean(abs(displ-spektraf));
    outkappa(2,i)=defr(i);
    end
    index5=find(abs(defr)==min(abs(defr)));
    kappa=kappa(index5);

%--------------------------------------------------------------------------
%                Calculation of Source Parameters
%--------------------------------------------------------------------------
m0=(4*pi*2.8*(3.5*100000)^3*omg)/(2*0.6);
mw=(0.6667*log10(m0))-10.7;
r=(2.34*(3.5*100000)/(2*pi*fc))/100;
sd=(7*m0)/(16*(r*100)^3)*0.000001;
m0=m0/10000000;%dyne-cm to Nm
sd=sd/10;%bars to MPa

%--------------------------------------------------------------------------
f=fv;lf=length(f);
% -------------------------------------------------------------------------
%                Plotting in acceleration spektra
% -------------------------------------------------------------------------
%               fmax
spektra=((2*pi*f).^2)*omg./((1.+(f./fc).^2).*(1.+(f./fmax).^p));% Brune Model with high cut
spektra1=((2*pi*fc).^2)*omg./((1.+(fc./fc).^2).*(1.+(fc./fmax).^p));
spektra2=((2*pi*fmax).^2)*omg./((1.+(fmax./fc).^2).*(1.+(fmax./fmax).^p));

%               kappa
spektrak=((2*pi*f).^2)*omg./((1.+(f./fc).^2).*exp(pi*fv*kappa));
spektrak1=((2*pi*fc).^2)*omg./((1.+(fc./fc).^2).*exp(pi*fc*kappa));
%               plotting for fmax

loglog(f,spektra,'b','LineWidth',2);grid on;
 plot(fc,spektra1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektra1,'r+', 'MarkerSize',25); 

 loglog(f,spektra,'b','LineWidth',2);grid on;
 plot(fc,spektra1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fmax,spektra2,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektra1,'r+', 'MarkerSize',25); plot(fmax,spektra2,'r+', 'MarkerSize',25);
 %              plotting for kappa
 loglog(f,spektrak,'k--','LineWidth',2);grid on;
 plot(fc,spektrak1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektrak1,'r+', 'MarkerSize',25); 

mostr=num2str(m0, '%3.1e\n');fcstrr=num2str(r, '%3.1f');
omstr=num2str(omg, '%3.1e\n');mstrp=num2str(p,'%3.1f');
mstr1=num2str(mw,'%3.1f');fcstr1=num2str(fc, '%3.1f');
mstr2=num2str(sd,'%5.3f');fcstr2=num2str(fmax, '%3.1f');
strk=num2str(kappa, '%4.3f');

tx=fv(1);
ty=min(spektra);

text(tx,ty*3,[' M_0= ',mostr,' Nm;','  {\it f}_c= ',fcstr1,' Hz;',...
    ' {\it f}_m_a_x= ', fcstr2,' Hz;',' p = ' ,mstrp,'; \kappa = ' ,strk],'HorizontalAlignment','left',...
    'BackgroundColor',[1 1 .9],'EdgeColor','r','LineWidth',1,'FontSize',10); 
text(tx,ty,[' M_w= ',mstr1,';','{\it r} = ',fcstrr,' m;', '\Delta\sigma=',mstr2,' MPa'],...
    'HorizontalAlignment','left','BackgroundColor',[1 1 .9],'EdgeColor','r',...
        'LineWidth',1,'FontSize',10);
    
saveas(gcf,['accspec' num2str(fff)],'fig');
saveas(gcf,['accspec' num2str(fff)],'bmp');hold off;close
% -------------------------------------------------------------------------
%                     PLOTTING OF DISPLACEMENT SPECTRA
% -------------------------------------------------------------------------
loglog(fv,displ,'m');grid on;title (['Displacement Spectrum:  ',num2str(fff)]);
xlabel('Frequency (Hz)'); ylabel('Displacement Spectrum');grid on;hold on
%               fmax
spektra=((2*pi*f).^0)*omg./((1.+(f./fc).^2).*(1.+(f./fmax).^p));% Brune Model with high cut
spektra1=((2*pi*fc).^0)*omg./((1.+(fc./fc).^2).*(1.+(fc./fmax).^p));
spektra2=((2*pi*fmax).^0)*omg./((1.+(fmax./fc).^2).*(1.+(fmax./fmax).^p));

%               kappa
spektrak=((2*pi*f).^0)*omg./((1.+(f./fc).^2).*exp(pi*fv*kappa));
spektrak1=((2*pi*fc).^0)*omg./((1.+(fc./fc).^2).*exp(pi*fc*kappa));
%               plotting for fmax

 loglog(f,spektra,'b','LineWidth',2);grid on;
 plot(fc,spektra1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektra1,'r+', 'MarkerSize',25); 

 loglog(f,spektra,'b','LineWidth',2);grid on;
 plot(fc,spektra1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fmax,spektra2,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektra1,'r+', 'MarkerSize',25); plot(fmax,spektra2,'r+', 'MarkerSize',25);
 %              plotting for kappa
 loglog(f,spektrak,'k--','LineWidth',2);grid on;
 plot(fc,spektrak1,'ro', 'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','g');
 plot(fc,spektrak1,'r+', 'MarkerSize',25); 

tx=fv(1);
ty=min(spektra);

text(tx,ty*3,[' M_0= ',mostr,' Nm;','  {\it f}_c= ',fcstr1,' Hz;',...
    ' {\it f}_m_a_x= ', fcstr2,' Hz;',' p = ' ,mstrp,'; \kappa = ' ,strk],'HorizontalAlignment','left',...
    'BackgroundColor',[1 1 .9],'EdgeColor','r','LineWidth',1,'FontSize',10); 
text(tx,ty,[' M_w= ',mstr1,';','{\it r} = ',fcstrr,' m;', '\Delta\sigma=',mstr2,' MPa'],...
    'HorizontalAlignment','left','BackgroundColor',[1 1 .9],'EdgeColor','r',...
        'LineWidth',1,'FontSize',10);
%--------------------------------------------------------------------------
%                        SAVING OF RESULTS
%--------------------------------------------------------------------------
saveas(gcf,['dispspec' num2str(fff)],'fig');
saveas(gcf,['dispspec' num2str(fff)],'bmp');close;
%--------------------------------------------------------------------------
load OUT.txt;
OUT(fff,1)=fff;
OUT(fff,2)=omg;
OUT(fff,3)=fc;
OUT(fff,4)=fmax;
OUT(fff,5)=p;
OUT(fff,6)=kappa;
OUT(fff,7)=m0;
OUT(fff,8)=mw;
OUT(fff,9)=r;
OUT(fff,10)=sd;
save OUT.txt  OUT* -ascii;   
clc;clear all
end
%Rows 1-9