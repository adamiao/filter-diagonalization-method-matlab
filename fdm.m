%-------------------------------------------------------------------------%
%
% This code was made based on the paper by Dr. Vladimir Mandelshtam:
% 'Harmonic inversion of time signals and its applications', Journal of
% Chemical Physics 107, 6756 (1997).
% 
% This code will, from a given 1D signal, decompose it into:
% 
%           ----
%           \
%    s(t) =  )    D_n * exp(-i (2 pi * F_n * t + P_n)) * exp(G_n * t)
%           /
%           ----
% 
% where D_n is the n-th amplitude content, F_n is the n-th frequency
% content, P_n is the n-th phase content and G_n is the n-th growth rate
% content.
%
%-------------------------------------------------------------------------%
%         FILTER DIAGONALIZATION METHOD - FOURIER-TYPE KRYLOV BASIS       %
%-------------------------------------------------------------------------%
%
% INPUT
% =====
% 'signal'    : 1-D signal (column or row)
% 'Fs'        : sampling frequency of the signal. It is equal to the number
%               of points of the signal divided by the total time elapsed
% 'fmin'      : minimum frequency being investigated (Hz)
% 'fmax'      : maximum frequency being investigated (Hz)
% 'eps'       : criteria for selecting the eigenvalues (1e-2 is suggested)
%
% SUGGESTION: start by using a small number of points (~300-500) as well as
%             a small frequency window (fmax-fmin)~5. Change these conditions
%             accordingly depending on how fast the code runs.
%
% OUTPUT
% ======
% 'freq'      : frequency of oscillation (Hz)
% 'growth'    : growth rate (positive means growth)
% 'amp'       : real amplitude
% 'phase'     : phase angle (degrees)
%
%-------------------------------------------------------------------------%

function [freq,growth,amp,phase]=fdm(signal,Fs,fmin,fmax,eps)


% CREATION OF MATRICES U and S
% ============================
% Given a value of f_min and f_max we have to stipulate the number of
% possible eigenvalues in this interval such that:
% 
%               2*pi*f_min*tau < phi_j < 2*pi*f_max*tau

L=length(signal);
tau=1/Fs;
wmin=2*pi*fmin;
wmax=2*pi*fmax;
M=floor((L-3)/2);
J=ceil(L*tau*(wmax-wmin)/(4*pi));
if J<4
    J=4;
end
disp('FDM decomposition');
disp(['   Freq Range = ' num2str(fmin) '-' num2str(fmax)]);
disp(['   Data Length L = ' num2str(L)]);
disp(['   Tau = 1/FS = ' num2str(tau)]);
disp(['   M = (L-3)/2 = ' num2str(M)]);
disp(['   J = L*dF/2FS = ' num2str(J)]);

phi=zeros(J,1);

delta=tau*(wmax-wmin)/(J-1);

for j=1:J
    phi(j)=tau*wmin+(j-1)*delta;
end

%disp('  FDM MATRIX');
[A,B,C]=fdm_matrix(signal,M,J,phi);

% GENERALIZED EIGENVALUE PROBLEM - AMPLITUDE AND FREQUENCY CALCULATION
% ====================================================================

%disp('  EIGEN VAL-VECT');
[E,D]=eig(B,A);

% % NORMALIZATION OF EIGENVECTORS
% Here we have to renormalize the vectors by means of the U^(0) matrix
% (here represented by A). This means that the renormalization norm is
% equal to: renorm_i = sqrt(eig_i^T * A * eig_i)

%disp('  NORMALIZE');
for i=1:J
    renorm=conj(E(:,i)')*A*E(:,i);
    E(:,i)=E(:,i)/sqrt(renorm);
end

% % FREQUENCY CALCULATION % %
%disp('  FREQUENCY CALCULATION');
f=[];
lambda=[];
V=[];
counter=0;
for i=1:J
    a=norm((C-D(i,i)^2*A)*E(:,i));
    b=1i/tau*log(D(i,i));
    if a<eps && real(b)>0 && real(b)<pi/tau
        counter=counter+1;
        f(counter,1)=b;
        lambda(counter,1)=D(i,i);
        V(:,counter)=E(:,i);
    end
end

% % AMPLITUDE CALCULATION % %
%disp('  AMPLITUDE CALCULATION');
d=fdm_amplitude(signal,M,J,phi,V,lambda,2);

[~,beta]=sort(abs(d),'descend');
d=d(beta);
f=f(beta);

% % ORGANIZATION OF OUTPUT DATA % %
%disp('  ORGANIZE OUTPUT DATA');
freq=real(f)/(2*pi);
growth=imag(f);
amp=abs(d);
phase=-atan2(imag(d),real(d))*180/pi;

disp(['   Num.EF = ' num2str(length(freq))]);

return

%-------------------------------------------------------------------------%
%          FILTER DIAGONALIZATION METHOD - AMPLITUDE CALCULATION          %
%-------------------------------------------------------------------------%
%
% INPUT
% =====
% 'signal'  : 1-D signal (column or row)
% 'M'       : number of points being considered in the signal
% 'J'       : size of matrices (generalized eigenvalue problem)
% 'phi'     : frequency grid that comes from FDM
% 'V'       : selected eigenvectors (out of the total J)
% 'lambda'  : selected eigenvalues (out of the total J)
% 'method'  : which method of calculation will be used (either 1 or 2). I
%             suggest method 2 be used: not much less accurate than method
%             1 and yet more reliable.
%
% OUTPUT
% ======
% 'd'       : complex amplitude
%
%-------------------------------------------------------------------------%
function [d]=fdm_amplitude(signal,M,J,phi,V,lambda,method)

% I suggest always using the second method. The first method has not been
% giving results as consistently as the second.

% % METHOD #1
% % =========
if method==1   
    [a,~]=size(lambda);
    d=zeros(a,1);
    for k=1:a
        for j=1:J
            U=fdm_amplitude_Uo(signal,M,exp(-1i*phi(j)),lambda(k));
            d(k)=d(k)+V(j,k)*U;
        end
        d(k)=(d(k)/(M+1))^2;
    end    
end

% % METHOD #2
% % =========
if method==2    
    [a,~]=size(lambda);
    d=zeros(a,1);
    for k=1:a
        for j=1:J            
            thePhiJ = 1i*phi(j);
            nvect = (0:M); 
            expThePhiJnVect = exp(thePhiJ*nvect);
            theFactor = dot(signal(1:length(expThePhiJnVect)),expThePhiJnVect);
            d(k)=d(k)+theFactor.*V(j,k);
        end
        d(k)=d(k)^2;
    end  
end

return

%-------------------------------------------------------------------------%
%                AMPLITUDE CALCULATION - AUXILIARY FUNCTION               %
%-------------------------------------------------------------------------%
% 
% Aux function used in conjunction with 'METHOD 1' of amplitude calculation.
% 
%-------------------------------------------------------------------------%

function [U]=fdm_amplitude_Uo(signal,M,z1,z2)

if abs(z1-z2)<1e-4 && abs(real(z1)-real(z2))<1e-4    
    sum=0;    
    for n=0:2*M
        sum=sum+(M-abs(M-n)+1)*signal(n+1)*z1;
    end   
    U=sum;    
else    
    U1=0; U2=0; U3=0; U4=0;  
    for n=0:M
        U1=signal(n+1)*z2^(-n)+U1;
    end
    U1=z1*U1;    
    for n=0:M
        U2=signal(n+1)*z1^(-n)+U2;
    end
    U2=-z2*U2;    
    for n=M+1:2*M
        U3=signal(n+1)*z2^(M-n+1)+U3;
    end
    U3=-z1^(-M)*U3;    
    for n=M+1:2*M
        U4=signal(n+1)*z1^(M-n+1)+U4;
    end
    U4=z2^(-M)*U4;    
    U=1/(z1-z2)*(U1+U2+U3+U4);   
end

return

%-------------------------------------------------------------------------%
%             FILTER DIAGONALIZATION METHOD - MATRIX CREATION             %
%-------------------------------------------------------------------------%
% This Matrix is symmetric (checked mathematically and numerically) 
%
% INPUT
% =====
% 'signal'  : 1-D signal (column or row)
% 'M'       : # terms in signal being used for correlation (<=L(signal)/2)
% 'J'       : size of the matrices
% 'phi'     : spatial grid where we want to look for particular frequencies
%
% OUTPUT
% ======
% 'A'       : is equivalent to the matrix U^(0)
% 'B'       : is equivalent to the matrix U^(1)
% 'C'       : is equivalent to the matrix U^(2)
%
%-------------------------------------------------------------------------%
function [A,B,C]=fdm_matrix(signal,M,J,phi)

A = zeros(J,J);
B = zeros(J,J);
C = zeros(J,J);

%disp('    not diagonal terms');
for i=1:J
    for j=i:J
        
        A1=0; A2=0; A3=0; A4=0;
        B1=0; B2=0; B3=0; B4=0;
        C1=0; C2=0; C3=0; C4=0;
        
        if i==j
            continue
        end
        
        thePhiI = 1i*phi(i);
        thePhiJ = 1i*phi(j);
        expThePhiI = exp(-thePhiI);
        expThePhiJ = exp(-thePhiJ);
        expThePhiIM = exp(thePhiI*M);
        expThePhiJM = exp(thePhiJ*M);
        invdifexpThePhiIJ = 1/(expThePhiI-expThePhiJ);
        
        nvect = (0:M);
        mvect = (M+1:2*M);  
        expThePhiJnVect = exp(thePhiJ*nvect);
        expThePhiInVect = exp(thePhiI*nvect);
        expThePhiJmVect = exp(-thePhiJ*(M-mvect+1));
        expThePhiImVect = exp(-thePhiI*(M-mvect+1));
        A1 = dot(signal(1:length(expThePhiJnVect)),expThePhiJnVect);
        B1 = dot(signal(2:length(expThePhiJnVect)+1),expThePhiJnVect);
        C1 = dot(signal(3:length(expThePhiJnVect)+2),expThePhiJnVect);
        A2 = dot(signal(1:length(expThePhiInVect)),expThePhiInVect);
        B2 = dot(signal(2:length(expThePhiInVect)+1),expThePhiInVect);
        C2 = dot(signal(3:length(expThePhiInVect)+2),expThePhiInVect);
        A3 = dot(signal(M+2:M+length(expThePhiJmVect)+1),expThePhiJmVect); 
        B3 = dot(signal(M+3:M+length(expThePhiJmVect)+2),expThePhiJmVect); 
        C3 = dot(signal(M+4:M+length(expThePhiJmVect)+3),expThePhiJmVect);
        A4 = dot(signal(M+2:M+length(expThePhiImVect)+1),expThePhiImVect); 
        B4 = dot(signal(M+3:M+length(expThePhiImVect)+2),expThePhiImVect); 
        C4 = dot(signal(M+4:M+length(expThePhiImVect)+3),expThePhiImVect);  
        A1=expThePhiI*A1;
        B1=expThePhiI*B1;
        C1=expThePhiI*C1;
        A2=-expThePhiJ*A2;
        B2=-expThePhiJ*B2;
        C2=-expThePhiJ*C2;    
        A3=-expThePhiIM*A3;
        B3=-expThePhiIM*B3;
        C3=-expThePhiIM*C3;     
        A4=expThePhiJM*A4;
        B4=expThePhiJM*B4;
        C4=expThePhiJM*C4;       
            
        A(i,j)=invdifexpThePhiIJ*(A1+A2+A3+A4);
        B(i,j)=invdifexpThePhiIJ*(B1+B2+B3+B4);
        C(i,j)=invdifexpThePhiIJ*(C1+C2+C3+C4);
        
        A(j,i)=A(i,j);
        B(j,i)=B(i,j);
        C(j,i)=C(i,j);  
     
    end
end

%disp('    diagonal terms');
for i=1:J      
    thePhiI = 1i*phi(i);
    nvect = (0:2*M); 
    factorexpThePhiInVect = (M-abs(M-nvect)+1).*exp(thePhiI*nvect);
    A(i,i) = dot(signal(1:length(factorexpThePhiInVect)),factorexpThePhiInVect);
    B(i,i) = dot(signal(2:length(factorexpThePhiInVect)+1),factorexpThePhiInVect);
    C(i,i) = dot(signal(3:length(factorexpThePhiInVect)+2),factorexpThePhiInVect);    
end

return
