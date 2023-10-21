function implizitschroedinger()
fprintf("Modus m=1: Elektron besitzt keine Geschwindigkeit und keine Wand. \n")
fprintf("Modus m=2: Elektron mit Geschwindigkeit und keiner Wand. \n")
fprintf("Modus m=3: Elektron mit Geschwindigkeit und Wand (Tunneleffekt) \n")
fprintf("Modus m=4: Elektron mit Geschwindigkeit und Wand mit Doppelspalt (Doppelspaltexperiment)\n")
L=input("Welcher Modus m soll gewaehlt werden? ");

n=30;

%Matrix A 
fprintf("Berechnet Matrix A...\n")
m=n^2;
f=zeros(m);
  A=zeros(m,m);
  for l=1:m
    for j=1:m
      if l == j
      A(l,j)=-2;
      endif
      if j-1==l && mod(l, n)!=0
      A(l, j)=1;
      endif
      if j+1==l && mod(l, n)!=1
      A(l, j)=1;
      endif
    endfor
  endfor

  
%Matrix B
fprintf("Berechnet Matrix B...\n")
v=n; %v=sqrt(m)
f=zeros(m);
  B=zeros(m,m);
  for l=1:m
    for j=1:m
      if l == j
      B(l,j)=-2;
      endif
      if j-v==l
      B(l, j)=1;
      endif
      if j+v==l
      B(l, j)=1;
      endif
    endfor
  endfor


%V Matrix fuer Potenzial  
  V=zeros(n,n);
  for l=1:n
    for j=1:n
      if j==1
      V(l, j)=10^6;
      endif
      if j==n
        V(l, j)=10^6;
      endif
      if l==1
        V(l, j)=10^6;
      endif
      if l==n
        V(l, j)=10^6;
      endif
    endfor
  endfor 
  
  if L==3
    %Wand (ohne Loch) --> Tunneleffekt
    V((n-round(n/4)):(n-round(n/4)+round(n/20)), :)=23000;
  elseif L==4
    %Wand mit Spalten --> Doppelspaltexperiment
    V((n-round(n/4)):(n-round(n/4)+round(n/20)), 1:round(0.3*n))=23000;
    V((n-round(n/4)):(n-round(n/4)+round(n/20)), round(0.4*n):round(0.6*n))=23000;
    V((n-round(n/4)):(n-round(n/4)+round(n/20)), round(0.7*n):n)=23000;
  endif
  
  %Umformung der Matrix V
  K=zeros(m, m);
  t=1;
  for j=1:n
    for i=1:n
      K(t, t)=V(j, i);
      t+=1;
    endfor
  endfor
  V=K;
  
%Abstand zwischen Stuetzpunkten und Zeitabstand
  h1=1/(n-1);
  h2=1/(10^5);

  x=linspace(0, 1, n).';

  y=linspace(0, 1, n).';
  
  g=0.1;
  
%Je nach gewaehltem m ist Geschwindigkeit 0 oder 25
if L==1
  k=0;
elseif L==2 || L==3 || L==4
  k=25;
endif
  
%initialisieren (Startzuweisung)
  u=[];
  for j=1:n
    for d=1:n
      u(d, j)=1/(sqrt(2*pi)*g)*exp(-0.5*((x(d)-0.5)^2+(y(j)-0.5)^2)/(g^2))*(cos(2*pi*k*x(d))+1i*sin(2*pi*k*x(d)));
    endfor
  endfor

  step=0;
  
  
%Matrix C
C=eye(n^2)-1i*A*h2/(h1^2)-1i*B*h2/(h1^2)+1i*h2*V;

fprintf("Inverse berechnen...\n")
vC=inv(C);
 
%Umwandlung von u in Vektor 
      p=zeros(m, 1);
      t=1;
      for j=1:n
          for i=1:n
            p(t)=u(j, i);
            t+=1;
          endfor
        endfor
      u=p; 
 
    %Beginn der Hauptschleife
    for q=0:h2:1000

      step+=1;
      
      %Naechster Schritt impliziteuler
      u=vC*u;
      
      
      %plotting
      if mod(step, 20)==0
        
        %u wieder in Matrixform bringen (als Matrix K)
        K=zeros(n,n);
        t=1;
        for j=1:n
            for i=1:n
              K(j, i)=u(t);
              t+=1;
            endfor
        endfor
        
        %Betrag der Matrix
        v=abs(K);
        
        %Faerben der Waende
        if L==3
          v((n-round(n/4)):(n-round(n/4)+round(n/20)), :)=0;
        elseif L==4
          v((n-round(n/4)):(n-round(n/4)+round(n/20)), 1:round(0.3*n))=0;
          v((n-round(n/4)):(n-round(n/4)+round(n/20)), round(0.4*n):round(0.6*n))=0;
          v((n-round(n/4)):(n-round(n/4)+round(n/20)), round(0.7*n):n)=0;
        endif
        
        colormap(hot)
        imagesc(x, y, v)
        set(gca, 'yDir', 'normal')
        axis equal
        colorbar()
        pause(0.01)
         
      endif
      
    
    endfor

endfunction












