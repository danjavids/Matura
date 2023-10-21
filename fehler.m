function fehler()

  
  n=30;
  werte1=[];
  werte2=[];
  werte3=[];
  
  
  %ist die Anzahl Schritte, bei welcher geplottet werden soll
  zaehler=200;
  
  achse=[];
  
%Matrizen für das explizite Euler-Verfahren
  
%Matrix A1   
  A1=zeros(n,n);
  for l=1:n
    for j=1:n
      if l == j
      A1(l,j)=-2;
      endif
      if j-1==l
      A1(l, j)=1;
      endif
      if j+1==l
      A1(l, j)=1;
      endif
    endfor
  endfor

%V1 Matrix fuer Potenzial  

  V1=zeros(n,n);
  for l=1:n
    for j=1:n
      if j==1
      V1(l, j)=10^6;
      endif
      if j==n
        V1(l, j)=10^6;
      endif
      if l==1
        V1(l, j)=10^6;
      endif
      if l==n
        V1(l, j)=10^6;
      endif
    endfor
  endfor 
  
  
%Matrizen für das implizite Euler-Verfahren

%Matrix A2 
fprintf("Berechnet Matrix A2...\n")
m=n^2;
f=zeros(m);
  A2=zeros(m,m);
  for l=1:m
    for j=1:m
      if l == j
      A2(l,j)=-2;
      endif
      if j-1==l && mod(l, n)!=0
      A2(l, j)=1;
      endif
      if j+1==l && mod(l, n)!=1
      A2(l, j)=1;
      endif
    endfor
  endfor

  
%Matrix B2
fprintf("Berechnet Matrix B2...\n")
v=n; %v=sqrt(m)
f=zeros(m);
  B2=zeros(m,m);
  for l=1:m
    for j=1:m
      if l == j
      B2(l,j)=-2;
      endif
      if j-v==l
      B2(l, j)=1;
      endif
      if j+v==l
      B2(l, j)=1;
      endif
    endfor
  endfor


%V2 Matrix fuer Potenzial  
  V2=zeros(n,n);
  for l=1:n
    for j=1:n
      if j==1
      V2(l, j)=10^6;
      endif
      if j==n
        V2(l, j)=10^6;
      endif
      if l==1
        V2(l, j)=10^6;
      endif
      if l==n
        V2(l, j)=10^6;
      endif
    endfor
  endfor
 
  %Umformung der Matrix V2
  K=zeros(m, m);
  t=1;
  for j=1:n
    for i=1:n
      K(t, t)=V2(j, i);
      t+=1;
    endfor
  endfor
  V2=K; 
  
%Abstand zwischen Stützpunkten und Zeitabstand
  h1=1/(n-1);
  h2=1/(5*10^4);
  
  x=linspace(0, 1, n).';

  y=linspace(0, 1, n).';
  
  g=0.1;
  
  k=0;
  
%initialisieren (Startzuweisung)
  u1=[];
  for j=1:n
    for d=1:n
      u1(d, j)=1/(sqrt(2*pi)*g)*exp(-0.5*((x(d)-0.5)^2+(y(j)-0.5)^2)/(g^2))*(cos(2*pi*k*x(d))+1i*sin(2*pi*k*x(d)));
    endfor
  endfor
  
  u2=u1;
  u3=u1;
  


  step=0;
  
%Matrix C
C=eye(n^2)-1i*A2*h2/(h1^2)-1i*B2*h2/(h1^2)+1i*h2*V2;

fprintf("Inverse berechnen...\n")
vC=inv(C);

%Umwandlung von u2 in Vektor 
      p=zeros(m, 1);
      t=1;
      for j=1:n
          for i=1:n
            p(t)=u2(j, i);
            t+=1;
          endfor
        endfor
      u2=p;   
 
    
    %Beginn der Hauptschleife
    for q=0:h2:1

      step+=1;
      
      %Schritt des expliziten Euler-Verfahrens
      
      u1=u1+1i*h2*(1/(h1*h1)*(u1*A1+A1*u1))-h2*1i*V1.*u1;
     
      u1(1, :)=0;
      u1(n, :)=0;
      u1(:, 1)=0;
      u1(:, n)=0;
      
      %Schritt des impliziten Euler-Verfahren
      u2=vC*u2;
      
      %Schritt des RK4
      k1=1i*h2*(1/(h1*h1)*(u3*A1+A1*u3))-h2*1i*V1.*u3;
      k2=1i*h2*(1/(h1*h1)*((u3+k1/2)*A1+A1*(u3+k1./2)))-h2*1i*V1.*(u3+k1/2);
      k3=1i*h2*(1/(h1*h1)*((u3+k2/2)*A1+A1*(u3+k2./2)))-h2*1i*V1.*(u3+k2/2);
      k4=1i*h2*(1/(h1*h1)*((u3+k3)*A1+A1*(u3+k3)))-h2*1i*V1.*(u3+k3);
      
      t=1/6*(k1+2*k2+2*k3+k4);
      
      u3=u3+t;
      
      u3(1, :)=0;
      u3(n, :)=0;
      u3(:, 1)=0;
      u3(:, n)=0;
      
      %plotting
      if mod(step, zaehler)==0
        %Betrag der Matrix u1
        v1=abs(u1);
        
       %u2 wieder in Matrixform bringen (als Matrix K)
        K=zeros(n,n);
        t=1;
        for j=1:n
            for i=1:n
              K(j, i)=u2(t);
              t+=1;
            endfor
        endfor
        
        
        
        
        %Betrag der Matrix u2/K
        v2=abs(K);
        
        %Betrag der Matrix u3
        v3=abs(u3);
        
        %Zusammenrechnen der Werte unter dem Graph (explizit)
        s1=0;
        for i=1:n
          for j=1:n
            s1=s1+v1(i,j).^2;
          endfor
        endfor
        s1=sqrt(s1);
        
        %Zusammenrechnen der Werte unter dem Graph (implizit)
        s2=0;
        for i=1:n
          for j=1:n
            s2=s2+v2(i,j).^2;
          endfor
        endfor
        s2=sqrt(s2);
        
        %Zusammenrechnen der Werte unter dem Graph (RK4)
        s3=0;
        for i=1:n
          for j=1:n
            s3=s3+v3(i,j).^2;
          endfor
        endfor
        s3=sqrt(s3);
        
        %Erster Wert als 100% definieren
        if step==zaehler
          p1=s1;
          p2=s2;
          p3=s3;
        endif
        prozent_explizit=s1/p1*100;
        prozent_implizit=s2/p2*100;
        prozent_RK4=s3/p3*100;
 
      %x-Achse: Zeitabstand
      %Prozentuale Abweichung vom Orginalwert 100%
      a=step/zaehler;
      achse(a)=step;
      werte1(a)=prozent_explizit;
      werte2(a)=prozent_implizit;
      werte3(a)=prozent_RK4;
      plot(achse, werte1, 'r')
      hold on
      plot(achse, werte2, 'b')
      hold on
      plot(achse, werte3, 'g')
      hold on 
      axis([0 step 0 300])
      grid on
      pause(0.001)
      
      if abs(s3-p3)>abs(s2-p2) && abs(s2-p2)>abs(s1-p1)
        fprintf('1.RK 2.implizit 3.explizit\n\n')
      endif
      if abs(s3-p3)>abs(s1-p1) && abs(s1-p1)>abs(s2-p2)
        fprintf('1.RK 2.explizit 3.implizit\n\n')
      endif
      if abs(s1-p1)>abs(s3-p3) && abs(s3-p3)>abs(s2-p2)
        fprintf('1.explizit 2.Rk 3.implizit\n\n')
      endif
      if abs(s1-p1)>abs(s2-p2) && abs(s2-p2)>abs(s3-p3)
        fprintf('1.explizit 2.implizit 3.Rk\n\n')
      endif
      if abs(s2-p2)>abs(s3-p3) && abs(s3-p3)>abs(s1-p1)
        fprintf('1.implizit 2.Rk 3.explizit\n\n')
      endif
      if abs(s2-p2)>abs(s1-p1) && abs(s1-p1)>abs(s3-p3)
        fprintf('1.implizit 2.explizit 3.RK\n\n')
      endif
      
    
      
      
      
     endif 
     if step==10000
        break
     endif
      
    endfor
  
  
endfunction