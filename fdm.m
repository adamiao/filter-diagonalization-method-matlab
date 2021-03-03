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
%
% 'signal'    : 1-D signal (column or row)
% 'Fs'        : sampling frequency of the signal. It is equal to the number
%               of points of the signal divided by the total time elapsed
% 'fmin'      : minimum frequency being investigated (Hz)
% 'fmax'      : maximum frequency being investigated (Hz)
% 'eps'       : criteria for selecting the eigenvalues (1e-2 is suggested)
%
% SUGGESTION: start by using a small number of points (~300-500) as well as a
%             small frequency window (fmax-fmin)~5. Change these conditions
%             accordingly depending on how fast the code runs.
%
% OUTPUT
% ======
%
% 'freq'      : frequency of oscillation (Hz)
% 'growth'    : growth rate (positive means growth)
% 'amp'       : real amplitude
% 'phase'     : phase angle (degrees)
%
%-------------------------------------------------------------------------%

function [freq,growth,amp,phase]=fdm(signal,Fs,fmin,fmax,eps)


% CREATION OF MATRICES U and S
% ============================
% 
% Given a value of f_min and f_max we have to stipulate the number of
% possible eigenvalues in this interval such that:
% 
%               2*pi*f_min*tau < phi_j < 2*pi*f_max*tau

L=length(signal); tau=1/Fs;
wmin=2*pi*fmin; wmax=2*pi*fmax;
M=floor((L-3)/2); J=ceil(L*tau*(wmax-wmin)/(4*pi));
if J<4
    J=4;
end

phi=zeros(J,1);

delta=tau*(wmax-wmin)/(J-1);

for j=1:J
    phi(j)=tau*wmin+(j-1)*delta;
end

[A,B,C]=fdm_matrix(signal,M,J,phi);

% GENERALIZED EIGENVALUE PROBLEM - AMPLITUDE AND FREQUENCY CALCULATION
% ====================================================================

[E,D]=eig(B,A);

% % NORMALIZATION OF EIGENVECTORS
% Here we have to renormalize the vectors by means of the U^(0) matrix
% (here represented by A). This means that the renormalization norm is
% equal to: renorm_i = sqrt(eig_i^T * A * eig_i)

for i=1:J
    renorm=conj(E(:,i)')*A*E(:,i);
    E(:,i)=E(:,i)/sqrt(renorm);
end

% % FREQUENCY CALCULATION % %

f=[]; lambda=[]; V=[]; counter=0;
for i=1:J
    a=norm((C-D(i,i)^2*A)*E(:,i));
    b=1i/tau*log(D(i,i));
    if a<eps && real(b)>0 && real(b)<pi/tau
        counter=counter+1;
        f(counter,1)=b;
        lambda(counter,1)=D(i,i); V(:,counter)=E(:,i);
    end
end

% % AMPLITUDE CALCULATION % %

d=fdm_amplitude(signal,M,J,phi,V,lambda,2);

[~,beta]=sort(abs(d),'descend');
d=d(beta);
f=f(beta);

% % ORGANIZATION OF OUTPUT DATA % %

freq=real(f)/(2*pi);
growth=imag(f);
amp=abs(d);
phase=-atan2(imag(d),real(d))*180/pi;

return

%-------------------------------------------------------------------------%
%          FILTER DIAGONALIZATION METHOD - AMPLITUDE CALCULATION          %
%-------------------------------------------------------------------------%
%
% INPUT
% =====
%
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
%
% 'd'       : complex amplitude
%
%-------------------------------------------------------------------------%

function [d]=fdm_amplitude(signal,M,J,phi,V,lambda,method)

% % METHOD #1
% % =========
% 
% I suggest always using the second method. The first method has not been
% giving results as consistently as the second.

if method==1
    
    [a,~]=size(lambda); d=zeros(a,1);
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
    
    [a,~]=size(lambda); d=zeros(a,1);
    for k=1:a
        for j=1:J
            for n=0:M
                d(k)=d(k)+signal(n+1)*exp(1i*phi(j)*n).*V(j,k);
            end
        end
        d(k)=d(k)^2;
    end
    
end

return

%-------------------------------------------------------------------------%
%                AMPLITUDE CALCULATION - AUXILIARY FUNCTION               %
%-------------------------------------------------------------------------%
% 
% Auxiliary function used in conjunction with 'METHOD 1' of the amplitude
% calculation.
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
%
% INPUT
% =====
%
% 'signal'  : 1-D signal (column or row)
% 'M'       : number of terms in the signal being used for the correlation
% 'J'       : size of the matrices
% 'phi'     : spatial grid where we want to look for particular frequencies
%
% OUTPUT
% ======
%
% 'A'       : is equivalent to the matrix U^(0)
% 'B'       : is equivalent to the matrix U^(1)
% 'C'       : is equivalent to the matrix U^(2)
%
%-------------------------------------------------------------------------%

function [A,B,C]=fdm_matrix(signal,M,J,phi)

for i=1:J
    for j=1:J
        
        A1=0; A2=0; A3=0; A4=0;
        B1=0; B2=0; B3=0; B4=0;
        C1=0; C2=0; C3=0; C4=0;
        
        if i==j
            continue
        end
        
        for n=0:M
            A1=signal(n+1)*exp(1i*phi(j)*n)+A1;
            B1=signal(n+2)*exp(1i*phi(j)*n)+B1;
            C1=signal(n+3)*exp(1i*phi(j)*n)+C1;
        end
        A1=exp(-1i*phi(i))*A1;
        B1=exp(-1i*phi(i))*B1;
        C1=exp(-1i*phi(i))*C1;
        
        for n=0:M
            A2=signal(n+1)*exp(1i*phi(i)*n)+A2;
            B2=signal(n+2)*exp(1i*phi(i)*n)+B2;
            C2=signal(n+3)*exp(1i*phi(i)*n)+C2;
        end
        A2=-exp(-1i*phi(j))*A2;
        B2=-exp(-1i*phi(j))*B2;
        C2=-exp(-1i*phi(j))*C2;
        
        for n=M+1:2*M
            A3=signal(n+1)*exp(-1i*phi(j)*(M-n+1))+A3;
            B3=signal(n+2)*exp(-1i*phi(j)*(M-n+1))+B3;
            C3=signal(n+3)*exp(-1i*phi(j)*(M-n+1))+C3;
        end
        A3=-exp(1i*phi(i)*M)*A3;
        B3=-exp(1i*phi(i)*M)*B3;
        C3=-exp(1i*phi(i)*M)*C3;
        
        for n=M+1:2*M
            A4=signal(n+1)*exp(-1i*phi(i)*(M-n+1))+A4;
            B4=signal(n+2)*exp(-1i*phi(i)*(M-n+1))+B4;
            C4=signal(n+3)*exp(-1i*phi(i)*(M-n+1))+C4;
        end
        A4=exp(1i*phi(j)*M)*A4;
        B4=exp(1i*phi(j)*M)*B4;
        C4=exp(1i*phi(j)*M)*C4;
        
        A(i,j)=1/(exp(-1i*phi(i))-exp(-1i*phi(j)))*(A1+A2+A3+A4);
        B(i,j)=1/(exp(-1i*phi(i))-exp(-1i*phi(j)))*(B1+B2+B3+B4);
        C(i,j)=1/(exp(-1i*phi(i))-exp(-1i*phi(j)))*(C1+C2+C3+C4);
        
    end
end

for i=1:J
    
    sumA=0; sumB=0; sumC=0;
    for n=0:2*M
        sumA=sumA+(M-abs(M-n)+1)*signal(n+1)*exp(1i*phi(i)*n);
        sumB=sumB+(M-abs(M-n)+1)*signal(n+2)*exp(1i*phi(i)*n);
        sumC=sumC+(M-abs(M-n)+1)*signal(n+3)*exp(1i*phi(i)*n);
    end
    A(i,i)=sumA;
    B(i,i)=sumB;
    C(i,i)=sumC;
    
end

return
