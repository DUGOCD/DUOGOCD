clear all
clc
tic;
XIS=[10,20,50,100];
DIM=[2,3];
DDWEI=[500];
 FUN9=['LSMOP1';'LSMOP3';'LSMOP5';'LSMOP8';'LSMOP9'];
% FUN9=['DTLZ1';'DTLZ2';'DTLZ3';'DTLZ4';'DTLZ5';'DTLZ6';'DTLZ7'];

IGD23=[];
HV23=[];
for dim=1:2
    IGD512=[];
    HV512=[];
    for fun9=1:5
        for ddwei=1:4
            IGD30=[];
            HV30=[];
            for sans=1:5
[IGD,HV1]=EIEAPT(DIM,XIS,FUN9,dim,fun9,ddwei);
% figure; 
% plot(1:Generations,F1daxiao,'b*');
IGD30=[IGD30,IGD];
HV30=[HV30,HV1];
fprintf('当前是 %d 个目标，第 %d 分组第 %d 个测试问题的第 %d 次run\n', DIM(dim), XIS(ddwei),fun9,sans);
      toc;
aaa=toc;
            end
              
            IGDfun9(ddwei,:)=IGD30;
            HVfun9(ddwei,:)=HV30;
      
        end
        IGD512=[IGD512;IGDfun9];
        HV512=[HV512;HVfun9];
    end
    IGD23=[IGD23;IGD512];
    HV23=[HV23;HV512];
end
filename1=sprintf('jiaoAblation%sIGD23512.xlsx','EIEAPT' );
filename2=sprintf('jiaoAblation%sHV23512.xlsx','EIEAPT'  );
xlswrite(filename1,IGD23);
xlswrite(filename2,HV23);
