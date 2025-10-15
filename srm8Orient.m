function [imgsSegmentd,nClasses,time,nEdges,szSegmsQs]=...
    srm8Orient(image,orientim,Qlevels,g,minSZ,W1)
tic%Statistical Region Merging on a square lattice with 8-connectivity
[Ix,Iy,Idiag45,Idiag135] = srmImGrad8(image,W1);%8-connectivity
Ix=max(abs(Ix),[],3);Iy=max(abs(Iy),[],3);Ix(:,end)=[];Iy(end,:)=[];
Idiag45=max(abs(Idiag45),[],3);Idiag135=max(abs(Idiag135),[],3);
Idiag45(:,end)=[];Idiag45(end,:)=[];
Idiag135(1,:)=[];Idiag135(:,end)=[];
[~,index]=sort(abs([Iy(:);Ix(:);Idiag45(:);Idiag135(:)]));
clear Ix Iy ;
nQ=numel(Qlevels);imgsSegmentd=cell(nQ,1);%maps=cell(nQ,1);
szIm=size(image);szIm1=szIm(1);szIm2=szIm(2);
n=szIm1*szIm2;imFinal=zeros(szIm);
A=reshape(1:n,szIm(1:2));treerank=zeros(szIm(1:2));
szSegms=ones(szIm(1:2));imSeg=image;
nEdges=numel(index);%Building pairs
idx2=reshape(A(:,1:end-1),[],1);idx1=reshape(A(1:end-1,:),[],1);
idIdiag45=reshape(A(1:end-1,1:end-1),[],1);
idIdiag135=reshape(A(2:end,1:end-1),[],1);
pairs1=[idx1;idx2;idIdiag45;idIdiag135];
pairs2=[idx1+1;idx2+szIm1;idIdiag45+szIm1+1;idIdiag135+szIm1-1];
clear idIdiag45 idIdiag135 idx2 idx1;
nClasses=zeros(1,nQ);logdelta=2*log(6*n);angBd=sin(pi/4);%angBd=sin(5*pi/18);

for Q=Qlevels
    iter=find(Q==Qlevels);  
    for i=1:nEdges
        C1=pairs1(index(i));C2=pairs2(index(i)); 
        diC1Angle=orientim(C1);diC2Angle=orientim(C2);
        C1Point=[1+rem(C1-1,szIm1),ceil(C1/szIm1)];
        C2Point=[1+rem(C2-1,szIm1),ceil(C2/szIm1)];
        C1ToC2Vec=C2Point-C1Point;
        if C1ToC2Vec(2)<0, C1ToC2Vec=-C1ToC2Vec;end
C1ToC2Angle=acos(C1ToC2Vec(1)/norm(C1ToC2Vec));
        dfDirectionC1C2=sin(abs(diC1Angle-diC2Angle));
        dfDirectionC1=sin(abs(diC1Angle-C1ToC2Angle));
        dfDirectionC2=sin(abs(diC2Angle-C1ToC2Angle));
dfAng=(dfDirectionC1C2<angBd)&(dfDirectionC1<angBd)&(dfDirectionC2<angBd);
% dfAng=(dfDirectionC1C2<angBd);
        while (A(C1)~=C1 ); C1=A(C1); end
        while (A(C2)~=C2 ); C2=A(C2); end
        dR=(imSeg(C1)-imSeg(C2))^2;
        dG=(imSeg(C1+n)-imSeg(C2+n))^2;dB=(imSeg(C1+2*n)-imSeg(C2+2*n))^2;        
        dev=g^2*logdelta*(1/szSegms(C1)+1/szSegms(C2))/(2*Q);   
        dev=1*1*dev;%dev=10, 500, 200, 100, 50*dev;
        predicat=( (dR<dev) & (dG<dev) & (dB<dev) );               
if (((C1~=C2)&predicat & dfAng) | szSegms(C1)<=minSZ | szSegms(C2)<=minSZ),
            if treerank(C1)>treerank(C2)%Find the new root for both regions
                A(C2) = C1; reg=C1;regNot=C2;
            elseif treerank(C1) < treerank(C2)
                A(C1) = C2; reg=C2;regNot=C1;
            elseif C1 ~= C2
                A(C2) = C1; reg=C1;regNot=C2;
                treerank(C1) = treerank(C1) + 1;
            end
            if C1~=C2% Merge regions
                nreg=szSegms(C1)+szSegms(C2);
            imSeg(reg)=(szSegms(C1)*imSeg(C1)+szSegms(C2)*imSeg(C2))/nreg;
imSeg(reg+n)=(szSegms(C1)*imSeg(C1+n)+szSegms(C2)*imSeg(C2+n))/nreg;
imSeg(reg+2*n)=(szSegms(C1)*imSeg(C1+2*n)+szSegms(C2)*imSeg(C2+2*n))/nreg;
                szSegms(reg)=nreg;nClasses(1,iter)=nClasses(1,iter)+1;
                szSegms(regNot)=0;
            end
        end
    end
    while 1
        map_ = A(A) ;
        if isequal(map_,A) ; break ; end
        A = map_ ;
    end   
    for i=1:3
        imFinal(:,:,i)=imSeg(A+(i-1)*n);
    end
    imgsSegmentd{iter}=imFinal;
    szSegmsOutput=szSegms(:);szSegmsOutput=nonzeros(szSegmsOutput);
    szSegmsQs{iter}=szSegmsOutput';
end
for t=nQ:-1:1,
    tSum=0;
    for j=1:t,
        tSum=tSum+nClasses(1,j);
    end
    nClasses(1,t)=tSum;
end
clear imSeg imFinal pairs1 pairs2;
nClasses=n-nClasses;time=toc%End of the function srm8