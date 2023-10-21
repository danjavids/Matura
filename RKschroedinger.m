function explizitschroedinger()
  fprintf("Modus m=1: Elektron besitzt keine Geschwindigkeit und keine Wand. \n")
  fprintf("Modus m=2: Elektron mit Geschwindigkeit und keiner Wand. \n")
  fprintf("Modus m=3: Elektron mit Geschwindigkeit und Wand (Tunneleffekt) \n")
  fprintf("Modus m=4: Elektron mit Geschwindigkeit und Wand mit Doppelspalt (Doppelspaltexperiment)\n")
  m=input("Welcher Modus m soll gewaehlt werden? ");

  n=101;
  
%Matrix A  
  A=zeros(n,n);
  for l=1:n
    for j=1:n
      if l == j
      A(l,j)=-2;
      endif
      if j-1==l
      A(l, j)=1;
      endif
      if j+1==l
      A(l, j)=1;
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
  
%Abstand zwischen Stützpunkten und Zeitabstand
  h1=1/(n-1);
  h2=1/(10^7);
  
  x=linspace(0, 1, n).';

  y=linspace(0, 1, n).';
  
  g=0.1;
  
%Je nach gewähltem m ist Geschwindigkeit 0 oder 25
  if m==1
    k=0;
  elseif m==2 || m==3 || m==4
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
   
  if m==3
    %Wand (ohne Loch) --> Tunneleffekt
    V((n-round(n/4)):(n-round(n/4)+round(n/20)), :)=23000;
  elseif m==4
    %Wand mit Loch --> Doppelspaltexperiment
    V((n-round(n/4)):(n-round(n/4)+round(n/20)), 1:round(0.3*n))=23000;
    V((n-round(n/4)):(n-round(n/4)+round(n/20)), round(0.4*n):round(0.6*n))=23000;
    V((n-round(n/4)):(n-round(n/4)+round(n/20)), round(0.7*n):n)=23000;
  endif
    
    %Beginn der Hauptschleife
    for q=0:h2:1

      step+=1;
      
      k1=1i*h2*(1/(h1*h1)*(u*A+A*u))-h2*1i*V.*u;
      k2=1i*h2*(1/(h1*h1)*((u+k1/2)*A+A*(u+k1./2)))-h2*1i*V.*(u+k1/2);
      k3=1i*h2*(1/(h1*h1)*((u+k2/2)*A+A*(u+k2./2)))-h2*1i*V.*(u+k2/2);
      k4=1i*h2*(1/(h1*h1)*((u+k3)*A+A*(u+k3)))-h2*1i*V.*(u+k3);
      
      t=1/6*(k1+2*k2+2*k3+k4);
      
      u=u+t;
     
      u(1, :)=0;
      u(n, :)=0;
      u(:, 1)=0;
      u(:, n)=0;
       
      %plotting
      if mod(step, 2000)==0
        %Betrag der Matrix
        v=abs(u);
        if m==3
          v((n-round(n/4)):(n-round(n/4)+round(n/20)), :)=0;
        elseif m==4
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












